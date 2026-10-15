// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

// The prediction route for -machine:device on made-up machines: which
// provider MS2 and the LSTM models (RT, CCS) ask for, and whether the
// choice warns. Printed so the test can pin the whole table.

#include <odia/DeviceRoute.h>

#include <iostream>
#include <string>

int main()
{
  ODIA::DeviceCaps mac;
  mac.coreml_provider = true;
  ODIA::DeviceCaps cuda;
  cuda.cuda_provider = cuda.nvidia_gpu = cuda.cudnn = true;
  cuda.cudnn_detail = "cuDNN 9.1.0";
  ODIA::DeviceCaps no_cudnn = cuda;
  no_cudnn.cudnn = false;
  no_cudnn.cudnn_detail = "no libcudnn.so on the library path";
  ODIA::DeviceCaps no_gpu;
  no_gpu.cuda_provider = true;
  ODIA::DeviceCaps cpu_ort;
  cpu_ort.nvidia_gpu = true;
  const std::size_t min = 200000;
  const struct { const char* name; std::string device; bool require; const ODIA::DeviceCaps* caps; std::size_t n; } cases[] = {
    {"mac auto small", "auto", false, &mac, 1000},
    {"mac auto large", "auto", false, &mac, 5000000},
    {"mac cpu large", "cpu", false, &mac, 5000000},
    {"mac coreml small", "coreml", false, &mac, 1000},
    {"mac auto require", "auto", true, &mac, 1000},
    {"cuda auto", "auto", false, &cuda, 1000},
    {"cuda no cudnn auto", "auto", false, &no_cudnn, 1000},
    {"cuda no cudnn require", "auto", true, &no_cudnn, 1000},
    {"cuda no cudnn forced", "cuda", false, &no_cudnn, 1000},
    {"cuda cpu", "cpu", false, &cuda, 1000},
    {"cuda provider no gpu", "auto", false, &no_gpu, 1000},
    {"cpu ort with gpu", "auto", false, &cpu_ort, 1000},
  };
  const auto provider = [](bool gpu, bool coreml) { return gpu ? (coreml ? "CoreML" : "CUDA") : "CPU"; };
  for (const auto& c : cases)
  {
    const auto r = ODIA::chooseDeviceRoute(c.device, c.require, *c.caps, c.n, min);
    std::cout << c.name << ": MS2 " << provider(r.ms2_gpu, r.ms2_coreml) << ", RT/CCS "
              << provider(r.lstm_gpu, r.lstm_coreml) << (r.warning.empty() ? "" : ", warns")
              << (r.error.empty() ? "" : ", refused") << "\n";
  }
  return 0;
}
