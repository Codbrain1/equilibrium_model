#pragma once
#include<pybind11/embed.h>
#include<filesystem>
#include<string>
#include<Python.h>

namespace pyth = pybind11;

class PythonWrapper {
	pyth::scoped_interpreter guard{};
	pyth::module_ python_module;
public:
	PythonWrapper()
	{
		try {
			std::filesystem::path p = "./scripts/visual";		//path класс для хранения путей
			std::string folder = p.parent_path().string();		//parent_path - возвращает родительскую директорию
			std::string filename = p.stem().string();		//stem - возвращает имя файла
			

			pyth::module_ sys = pyth::module_::import("sys"); // добавляем путь в python sys.path
			sys.attr("path").attr("append")(folder);

			python_module = pyth::module_::import(filename.c_str());	//загружаем файл python
		}
		catch (const pyth::error_already_set& e) {
			throw std::runtime_error("Ошибка загрузки Python: " + std::string(e.what()));
		}
	}
	void vcl()
	{
		try {
			python_module.attr("visual_converation_laws")();
		}
		catch (const pyth::error_already_set& e) {
			throw std::runtime_error("Ошибка загрузки функции: " + std::string(e.what()));
		}
	}
	void vt()
	{
		try {
			python_module.attr("visual_traectories")();
		}
		catch(const pyth::error_already_set& e){
			throw std::runtime_error("Ошибка загрузки функции:" + std::string(e.what()));
		}
	}
	void ptc()
	{
		try {
			python_module.attr("print_to_cadrs")();
		}
		catch (const pyth::error_already_set& e) {
			throw std::runtime_error("Ошибка загрузки функции:" + std::string(e.what()));
		}
	}

};