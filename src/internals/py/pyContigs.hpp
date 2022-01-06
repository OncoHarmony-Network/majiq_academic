/**
 * pyContigs.hpp
 *
 * Copyright 2020 <University of Pennsylvania>
 */

#ifndef MAJIQ_PYBIND_PYCONTIGS_HPP
#define MAJIQ_PYBIND_PYCONTIGS_HPP

#include <pybind11/pybind11.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "../Contigs.hpp"

#include "templateTypes.hpp"

namespace majiq {
namespace bindings {

using pyContigs_t = pyClassShared_t<majiq::Contigs>;

inline void init_Contigs(pyContigs_t& pyContigs) {
  pyContigs
    .def("checksum",
        [](const Contigs& self) { return majiq::checksum(self); },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "checksum of contigs")
    .def(
        pybind11::init([](pybind11::list seqids) {
          auto result = Contigs::create();
          for (auto seqid_py : seqids) {
            auto seqid = seqid_py.cast<seqid_t>();
            if (result->count(seqid)) {
              throw std::invalid_argument("Contigs require unique seqids");
            }
            result->add(Contig{seqid});
          }
          return result;
        }),
        "Set up Contigs object using specified identifiers",
        pybind11::arg("seqids"))
    .def_property_readonly("seqid", &Contigs::seqids,
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        R"pbdoc(
        Sequence[str] of contig ids in order matching contig_idx
        )pbdoc")
    .def("__len__", &Contigs::size,
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def("__contains__",
        [](const Contigs& s, seqid_t x) -> bool { return s.count(x) > 0; },
        pybind11::call_guard<pybind11::gil_scoped_release>())
    .def("__getitem__",
        [](const Contigs& self, seqid_t x) { return self.get_idx(x); },
        pybind11::call_guard<pybind11::gil_scoped_release>(),
        "contig_idx for specified seqid", pybind11::arg("seqid"));
}

}  // namespace bindings
}  // namespace majiq

#endif  // MAJIQ_PYBIND_PYCONTIGS_HPP
