cmake_minimum_required(VERSION 3.10)
project(n-body LANGUAGES CXX)

# Настройка стандарта C++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
# Включение OpenMP (обязательно ДО add_executable)
find_package(OpenMP REQUIRED)  # Изменили на REQUIRED

# Список исходных файлов
set(SOURCES
    main.cpp
    functions.cpp
)

# Создание исполняемого файла
add_executable(n-body ${SOURCES})  # Используем имя n-body как в ошибке

# Подключение заголовков
target_include_directories(n-body PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Настройки компиляции
target_compile_options(n-body PRIVATE
    -march=native
    -O2
    -g
)

# Подключение OpenMP (обязательно ПОСЛЕ add_executable)
target_link_libraries(n-body PRIVATE OpenMP::OpenMP_CXX)

# Для clang может потребоваться явно указать библиотеку атомарных операций
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_libraries(n-body PRIVATE atomic)
endif()