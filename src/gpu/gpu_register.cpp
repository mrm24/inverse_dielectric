/// Copyright (C) 2020-2024 GreenX library
///
/// Licensed under the Apache License, Version 2.0 (the "License");
/// you may not use this file except in compliance with the License.
/// You may obtain a copy of the License at
///
///   http://www.apache.org/licenses/LICENSE-2.0
///
/// Unless required by applicable law or agreed to in writing, software
/// distributed under the License is distributed on an "AS IS" BASIS,
/// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
/// implied. See the License for the specific language governing
/// permissions and limitations under the License.
///
///  @file
///  This file contains a host-device register, which can in principle can be used to
///  create a fine map host-device map address (even for derived types). 
///  I have exposed C functions for Fortran interfacing.
///  Limits: each host must have its own register, though it can connect to several devices.
///  Dependencies: OMP4.5 or higher
///
#include <unordered_map>
#include <tuple>
#include <string>
#include <stdexcept>
#include <iostream>
#include <omp.h>

/// Alias for: {Device id, size, host_ptr, device_ptr}
using ref_info      = std::tuple<std::size_t, std::size_t, void*, void*>;
using identifier    = std::string;
using register_type = std::unordered_map<identifier, ref_info>;

class device_host_register {
private:
    /// Internal register
    register_type rg;
    /// Host id, this is obtained through OMP call
    const int host_id;
public:
    /// Delete the copy constructor.
    device_host_register(const device_host_register& source) = delete;
    /// Delete assignment operator.
    device_host_register& operator=(const device_host_register& source) = delete;

    /// Inits the register (the host identity is obtained through an omp call)
    device_host_register() :  host_id(omp_get_initial_device()) {}

    /// When out of scope or explicitly deallocated the pointer in C compliant calls
    /// It checks the registered addresses and cleans them
    ~device_host_register() {
        for (auto &entry : this->rg) {
            if (std::get<2>(entry.second) == nullptr) {
                omp_target_free(std::get<3>(entry.second), std::get<0>(entry.second));
            } else {
                auto info = omp_target_disassociate_ptr(std::get<2>(entry.second), std::get<0>(entry.second));
                omp_target_free(std::get<3>(entry.second), std::get<0>(entry.second));
            }
        }
    }

    /// Allocs memory on the device
    /// @param[in] id        - name identifying the variable
    /// @param[in] size      - the bytes that should be allocated in the device
    /// @param[in] device_id - device id in which the memory should be allocated
    void alloc_device(identifier id, std::size_t size, std::size_t device_id) {
        if (rg.count(id) != 0) {
            throw std::runtime_error("Error in alloc_device: " + id + " is already in use");
        }

        auto d_ptr = omp_target_alloc(size, device_id);

        if (d_ptr == nullptr) {
            throw std::runtime_error("Error in alloc: " + id + ". Call to omp_target_alloc failed");
        }

        rg[id] = {device_id, size, nullptr, d_ptr};
    }

    /// Associates a host pointer to a device pointer
    /// @param[in] id        - name identifying the variable
    /// @param[in] host_ptr  - the pointer in the host which should be associated to the id
    void associate(identifier id, void* host_ptr) {
        // Check existence
        if (rg.count(id) == 0) {
            throw std::runtime_error("Error in associate: " + id + " is not in the list");
        }
        auto& mrg = rg[id];
        // Check if the device pointer is already associated
        if (std::get<2>(mrg) != nullptr) {
            throw std::runtime_error("Error in associate: " + id + " device pointer is already associated to a host pointer");
        }
        // Check if the host pointer is associated
        if (omp_target_is_present(host_ptr, std::get<0>(mrg))) {
            throw std::runtime_error("Error in associate: " + id + " host pointer is already associated to a device pointer");
        }

        auto info = omp_target_associate_ptr(host_ptr, std::get<3>(mrg), std::get<1>(mrg), 0, std::get<0>(mrg));

        if (info != 0) {
            throw std::runtime_error("Error in associate: " + id + ". Call to omp_target_associate_ptr failed");
        }

        std::get<2>(mrg) = host_ptr;
    }

    /// Associates a host pointer to a device pointer
    /// @param[in] id        - name identifying the variable which should be disassociated
    void disassociate(identifier id) {
        // Check existence
        if (rg.count(id) == 0) {
            throw std::runtime_error("Error in disassociate: " + id + " is not in the list");
        }
        auto& mrg = rg[id];
        // Check if the device pointer is already associated
        if (std::get<2>(mrg) == nullptr) {
            throw std::runtime_error("Error in disassociate: " + id + " device pointer is not associated to a host pointer");
        }
        auto info = omp_target_disassociate_ptr(std::get<2>(mrg), std::get<0>(mrg));
        
        if (info != 0) {
            throw std::runtime_error("Error in disassociate: " + id + ". Call to omp_target_disassociate_ptr failed");
        }

        std::get<2>(mrg) = nullptr;
    }

    /// Send data between the host and the device, namely between memory addresses associated to the id
    /// @param[in] id        - name identifying the variable that should be transferred from the host to the device
    void host_to_device(identifier id) {
        // Check existence
        if (rg.count(id) == 0) {
            throw std::runtime_error("Error in host_to_device: " + id + " is not in the list");
        }
        auto& mrg = rg[id];
        // Check if everything is associated
        if (std::get<2>(mrg) == nullptr || std::get<3>(mrg) == nullptr) {
            throw std::runtime_error("Error in host_to_device: " + id + " some pointers are not associated");
        }

        auto info = omp_target_memcpy(std::get<3>(mrg), std::get<2>(mrg), std::get<1>(mrg), 0, 0, std::get<0>(mrg), this->host_id);

        if (info != 0) {
            throw std::runtime_error("Error in host_to_device: " + id + ". Call to omp_target_memcpy failed");
        }
    }

    /// Send data between the device and the host, namely Between memory addresses associated to the id
    /// @param[in] id  - name identifying the variable that should be transferred from the device to the host
    void device_to_host(identifier id) {
        // Check existence
        if (rg.count(id) == 0) {
            throw std::runtime_error("Error in device_to_host: " + id + " is not in the list");
        }
        auto& mrg = rg[id];
        // Check if everything is associated
        if (std::get<2>(mrg) == nullptr || std::get<3>(mrg) == nullptr) {
            throw std::runtime_error("Error in device_to_host: " + id + " some pointers are not associated");
        }
        
        auto info = omp_target_memcpy(std::get<2>(mrg), std::get<3>(mrg), std::get<1>(mrg), 0, 0, this->host_id, std::get<0>(mrg));

        if (info != 0) {
            throw std::runtime_error("Error in device_to_host: " + id + ". Call to omp_target_memcpy failed");
        }
    }

    /// Removes id from the register, it cleans and disassociate if required
    /// @param[in] id        - name identifying the variable that should be transferred from the host to the device
    void remove(identifier id) {
        if (rg.count(id) == 0) {
            throw std::runtime_error("Error in remove: " + id + " is not in the list");
        }
        auto& mrg = rg[id];
        if (std::get<2>(mrg) == nullptr) {
            omp_target_free(std::get<3>(mrg), std::get<0>(mrg));
        } else {
            auto info = omp_target_disassociate_ptr(std::get<2>(mrg), std::get<0>(mrg));
            if (info != 0) {
                throw std::runtime_error("Error in remove: " + id + ". Call to omp_target_disassociate_ptr failed");
            }
            omp_target_free(std::get<3>(mrg), std::get<0>(mrg));
        }
        rg.erase(id);
    }

    /// Get device pointer associated to an entry:
    /// @param[in] id        - name identifying the variable for which the device pointer will be retrieved
    void* get_device_ptr(identifier id) {
        if (rg.count(id) == 0) {
            throw std::runtime_error("Error in get_device_ptr: " + id + " is not in the list");
        }
        auto &mrg = rg[id];
        return std::get<3>(mrg);
    }

};

/// Expose C++ class and members for Fortran interfacing
extern "C" {

    void _constructor_device_host_register(device_host_register** rg) {
        try {
            *rg = new device_host_register;
        } catch (const std::bad_alloc& e) {
            throw std::runtime_error("Error in create_register : allocation of the register object failed");
        }
    }

    void _destructor_device_host_register(device_host_register* rg) {
        delete rg;
    }

    void _alloc_device_device_host_register(device_host_register* rg, const char* id, std::size_t size, int device_id) {
        std::string cxx_id(id);
        rg->alloc_device(cxx_id, size, device_id);
    }

    void _associate_device_device_host_register(device_host_register* rg, const char* id, void* host_ptr) {
        std::string cxx_id(id);
        rg->associate(cxx_id, host_ptr);
    }

    void _disassociate_device_device_host_register(device_host_register* rg, const char* id) {
        std::string cxx_id(id);
        rg->disassociate(cxx_id);
    }

    void _host_to_device_device_device_host_register(device_host_register* rg, const char* id){
        std::string cxx_id(id);
        rg->host_to_device(cxx_id);
    }

    void _device_to_host_device_device_host_register(device_host_register* rg, const char* id){
        std::string cxx_id(id);
        rg->device_to_host(cxx_id);
    }

    void _remove_device_device_host_register(device_host_register* rg, const char* id){
        std::string cxx_id(id);
        rg->remove(cxx_id);
    }
    void*  _get_device_ptr_device_device_host_register(device_host_register* rg, const char* id){
        std::string cxx_id(id);
        auto d_ptr = rg->get_device_ptr(cxx_id);
        return d_ptr;
    }

}

