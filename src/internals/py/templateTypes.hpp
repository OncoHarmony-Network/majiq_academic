/**
 * templateTypes.hpp
 *
 * python bindings to C++ classes using shared_ptr
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_PYBIND_TEMPLATETYPES_HPP
#define MAJIQ_PYBIND_TEMPLATETYPES_HPP

#include <pybind11/pybind11.h>

#include <memory>

namespace majiq {
namespace bindings {

template <typename T>
using pyClassShared_t = pybind11::class_<T, std::shared_ptr<T>>;

template <typename T>
using base_t = std::remove_cv_t<std::remove_reference_t<T>>;

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_TEMPLATETYPES_HPP
