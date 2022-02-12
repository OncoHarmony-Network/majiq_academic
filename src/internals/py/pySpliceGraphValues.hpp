/**
 * pySpliceGraphValues.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYSPLICEGRAPHVALUES_HPP
#define MAJIQ_PYBIND_PYSPLICEGRAPHVALUES_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <memory>
#include <utility>
#include <vector>

#include "../SpliceGraphReads.hpp"
#include "../SpliceGraphMask.hpp"

#include "templateTypes.hpp"
#include "vectorArray.hpp"

namespace majiq {
namespace bindings {

using pySpliceGraphReads_t = pyClassShared_t<SpliceGraphReads>;
using pySpliceGraphMask_t = pyClassShared_t<SpliceGraphMask>;

// we have expose_type, impl_type because std::vector<bool> is a non-container
// specialization because of C++ standards weirdness (packing 8 bools per
// byte). This breaks array views. So vector<impl_type> internally, but we want
// to consider its values as expose_type.
template <
  typename SGValuesT,
  typename expose_type = typename SGValuesT::value_type,
  typename impl_type = typename SGValuesT::value_type,
  std::enable_if_t<sizeof(expose_type) == sizeof(impl_type), bool> = true>
inline void define_splicegraphvalues_properties(
    pyClassShared_t<SGValuesT>& pySGValues) {
  pySGValues
    .def_property_readonly("_introns",
        &SGValuesT::introns,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Underlying introns")
    .def_property_readonly("_junctions",
        &SGValuesT::junctions,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Underlying junctions")
    .def_property_readonly("introns_values",
        [](pybind11::object& self_obj) {
        SGValuesT& self = self_obj.cast<SGValuesT&>();
        return ArrayFromVectorAndOffset<expose_type, impl_type>(
            self.introns_values(), 0, self_obj);
        },
        "values for each intron")
    .def_property_readonly("junctions_values",
        [](pybind11::object& self_obj) {
        SGValuesT& self = self_obj.cast<SGValuesT&>();
        return ArrayFromVectorAndOffset<expose_type, impl_type>(
            self.junctions_values(), 0, self_obj);
        },
        "values for each junction")
    .def(pybind11::init([](
            const std::shared_ptr<GeneIntrons>& introns,
            const std::shared_ptr<GeneJunctions>& junctions,
            pybind11::array_t<expose_type> _introns_values,
            pybind11::array_t<expose_type> _junctions_values) {
          if (_introns_values.ndim() != 1 || _junctions_values.ndim() != 1) {
            throw std::invalid_argument(
                "introns/junctions values must both be 1D");
          }
          std::vector<impl_type> ivalues_vec(_introns_values.shape(0));
          {
            auto introns_values = _introns_values.template unchecked<1>();
            for (size_t i = 0; i < ivalues_vec.size(); ++i) {
            ivalues_vec[i] = static_cast<impl_type>(introns_values(i));
            }
          }
          std::vector<impl_type> jvalues_vec(_junctions_values.shape(0));
          {
            auto junctions_values = _junctions_values.template unchecked<1>();
            for (size_t j = 0; j < jvalues_vec.size(); ++j) {
            jvalues_vec[j] = static_cast<impl_type>(junctions_values(j));
            }
          }
          pybind11::gil_scoped_release release;  // release GIL at this stage
          return SGValuesT{introns, junctions,
              std::move(ivalues_vec), std::move(jvalues_vec)};
          }),
        "Initialize values over splicegraph from numpy arrays",
        pybind11::arg("introns"),
        pybind11::arg("junctions"),
        pybind11::arg("introns_values"),
        pybind11::arg("junctions_values"));
}

inline void init_SpliceGraphReads(pySpliceGraphReads_t& pySpliceGraphReads) {
  define_splicegraphvalues_properties(pySpliceGraphReads);
  pySpliceGraphReads
    .def_static("from_sj", &SpliceGraphReads::FromSJ,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "Obtain raw readrates for introns/junctions from experiment SJ",
        pybind11::arg("introns"),
        pybind11::arg("junctions"),
        pybind11::arg("sj_introns"),
        pybind11::arg("sj_junctions"));
}

inline void init_SpliceGraphMask(pySpliceGraphMask_t& pySpliceGraphMask) {
  // explicitly set exposed type as bool (instead of underlying char)
  define_splicegraphvalues_properties<SpliceGraphMask, bool>(pySpliceGraphMask);
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYSPLICEGRAPHVALUES_HPP
