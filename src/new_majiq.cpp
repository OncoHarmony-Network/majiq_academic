/**
 * new_majiq.cpp
 *
 * Python bindings to new MAJIQ code
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void init_SpliceGraphAll(py::module_&);

PYBIND11_MODULE(internals, m) {
  m.attr("__name__") = "new_majiq.internals";
  // documentation of module
  m.doc() = R"pbdoc(
  MAJIQ classes for splicegraph

  .. current_module::new_majiq
  .. autosummary::
    :toctree: _generate

  )pbdoc";


  // vector<string> will be bound to Python list
  py::bind_vector<std::vector<std::string>>(m, "VectorString");

  init_SpliceGraphAll(m);


#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}  // PYBIND11_MODULE(new_majiq, m)
