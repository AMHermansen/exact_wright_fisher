#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <deque>
#include <sstream>
#include <stdexcept>
#include <vector>


#include "WrightFisher.h"

namespace py = pybind11;

double dt_default = 0.1;
double bt_default = 0.04;


double wf_end(WrightFisher& wf, double100 t){
     const Options o(dt_default, bt_default);
     return wf.DrawEndpoint(0., 0., t, o, wf.WF_gen).first;
}


vector<double> wf_end_points(WrightFisher& wf, double100 t, const int n_sims) {
     std::vector<double> wf_ends(n_sims);
     const Options o(dt_default, bt_default);
     for (int i = 0; i < n_sims; i++) {
          wf_ends[i] = wf_end(wf, t);
     }
     std::cout << wf_ends.size() << '\n';
     return wf_ends;
}


PYBIND11_MODULE(EWF_pybind, m) {
  py::class_<WrightFisher>(m, "WrightFisher")
      .def(py::init<vector<double>, bool, double100, int, double, int,
                    vector<double>>(),
           py::arg("thetaP"), py::arg("non_neut"), py::arg("sigma"),
           py::arg("selectionSetup"), py::arg("dom"), py::arg("SelPolyDeg"),
           py::arg("selCoefs"))
      .def("DiffusionRunner", &WrightFisher::DiffusionRunner, py::arg("nSim"),
           py::arg("x"), py::arg("startT"), py::arg("endT"),
           py::arg("Absorption"), py::arg("Filename"),
           py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default)
      .def("BridgeDiffusionRunner", &WrightFisher::BridgeDiffusionRunner,
           py::arg("nSim"), py::arg("x"), py::arg("z"), py::arg("startT"),
           py::arg("endT"), py::arg("sampleT"), py::arg("Absorption"),
           py::arg("Filename"), py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default)
      .def("DiffusionDensityCalculator",
           &WrightFisher::DiffusionDensityCalculator, py::arg("meshSize"),
           py::arg("x"), py::arg("startT"), py::arg("endT"),
           py::arg("Absorption"), py::arg("Filename"),
           py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default)
      .def("BridgeDiffusionDensityCalculator",
           &WrightFisher::BridgeDiffusionDensityCalculator, py::arg("meshSize"),
           py::arg("x"), py::arg("z"), py::arg("startT"), py::arg("endT"),
           py::arg("sampleT"), py::arg("Absorption"), py::arg("Filename"),
           py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default)
      .def("wright_fisher_sample", &wf_end, py::return_value_policy::copy)
      .def("wright_fisher_samples", [](WrightFisher& wf, double100 t, py::array_t<double> wf_end_points_array) {

           auto n_samples = wf_end_points_array.size();
           py::buffer_info buf = wf_end_points_array.request();
           auto ptr_buf = static_cast<double *>(buf.ptr);

           for (int i = 0; i < n_samples; i++) {
                ptr_buf[i] = wf_end(wf, t);
                std::cout << ptr_buf[i] << '\n';
           }
           return wf_end_points_array;
      });
}