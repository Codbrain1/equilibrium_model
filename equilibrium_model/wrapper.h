#pragma once
#include <pybind11/embed.h>

namespace py = pybind11;

class PythonWrapper {
public:
    PythonWrapper();
    ~PythonWrapper()=default;

    void visual_converation_laws();
    void visual_traectories();
    void vs();
    void vectors_velocity();
    void runPythonScript(const std::string& code);

private:
    py::scoped_interpreter guard{}; // Интерпретатор должен жить дольше всех вызовов
};