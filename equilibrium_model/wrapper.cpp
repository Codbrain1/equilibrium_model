#include "wrapper.h"
#include <stdexcept>
#include<iostream>
PythonWrapper::PythonWrapper() {
    // Добавляем путь к скриптам
    py::module_ sys = py::module_::import("sys");
    sys.attr("path").attr("append")("./scripts");
}

void PythonWrapper::visual_converation_laws() {
    try {
        std::cout << "Attempting to import visual module..." << std::endl;

        // Проверка существования файла
        auto os = py::module_::import("os");
        std::string script_path = "./scripts/visual.py";
        bool exists = os.attr("path").attr("exists")(script_path).cast<bool>();

        if (!exists) {
            throw std::runtime_error("File not found: " + script_path);
        }

        std::cout << "Found script at: " << script_path << std::endl;

        auto module = py::module_::import("visual"); // Без .py!
        std::cout << "Module imported successfully" << std::endl;

        module.attr("visual_converation_laws")();
    }
    catch (const py::error_already_set& e) {
        std::cerr << "Python error details:\n" << e.what() << std::endl;
        throw std::runtime_error("Python operation failed");
    }
    catch (const std::exception& e) {
        std::cerr << "Standard error: " << e.what() << std::endl;
        throw;
    }
}

void PythonWrapper::visual_traectories() {
    try {
        auto module = py::module_::import("visual");
        module.attr("visual_traectories")();
    }
    catch (const py::error_already_set& e) {
        throw std::runtime_error(e.what());
    }
}
void PythonWrapper::vs() {
    try {
        auto module = py::module_::import("visual");
        module.attr("vs")();
    }
    catch (const py::error_already_set& e) {
        throw std::runtime_error(e.what());
    }
}
void PythonWrapper::vectors_velocity() {
    try {
        auto module = py::module_::import("visual");
        module.attr("vectors_velocity")();
    }
    catch (const py::error_already_set& e) {
        throw std::runtime_error(e.what());
    }
}

void PythonWrapper::runPythonScript(const std::string& code) {
    py::exec(code);
}