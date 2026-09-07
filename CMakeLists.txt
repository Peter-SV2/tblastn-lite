# tblastn-lite.
#
# ImGui is FETCHED at configure time, pinned to a tag, rather than vendored:
# pinning is the point, and a tag gives it without putting fifty thousand
# lines of someone else's source in the history.
#
#   cmake -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build
#
# The search itself is a header with no window and no GPU in it, so the
# self-check builds and runs anywhere:
#
#   cmake --build build --target check
cmake_minimum_required(VERSION 3.16)
project(tblastn_lite CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

if(MSVC)
  add_compile_options(/W4 /permissive- /utf-8 $<$<CONFIG:Release>:/O2>)
else()
  add_compile_options(-Wall -Wextra)
endif()

if(MINGW)
  # Static for every target. Linked dynamically, the exe needs libstdc++-6.dll
  # beside it -- and if some other libstdc++ is earlier on PATH it loads that
  # one and faults inside <fstream> before main does anything. That is exactly
  # how the first release failed on the machine it was built on.
  add_link_options(-static -static-libgcc -static-libstdc++)
endif()

find_package(Threads REQUIRED)

add_executable(tblastn_lite tblastn_lite.cpp)
target_link_libraries(tblastn_lite PRIVATE Threads::Threads)

enable_testing()
add_test(NAME selftest COMMAND tblastn_lite --selftest)
add_custom_target(check COMMAND tblastn_lite --selftest DEPENDS tblastn_lite)

# --- the GUI ------------------------------------------------------------------
# Windows only: the host loop is Win32 + D3D11. Everything above is portable.
if(WIN32)
  include(FetchContent)
  FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.91.5
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(imgui)

  add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/misc/cpp/imgui_stdlib.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_win32.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_dx11.cpp)
  target_include_directories(imgui PUBLIC
    ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends ${imgui_SOURCE_DIR}/misc/cpp)
  # Third-party source is not ours to keep warning-clean.
  if(MSVC)
    target_compile_options(imgui PRIVATE /W0)
  else()
    target_compile_options(imgui PRIVATE -w)
  endif()

  add_executable(tblastn_gui WIN32 tblastn_imgui.cpp)
  # d3dcompiler because this backend compiles its shaders at runtime; without it
  # the link fails on D3DCompile and nothing says which library is missing.
  target_link_libraries(tblastn_gui PRIVATE imgui Threads::Threads d3d11 dxgi d3dcompiler comdlg32)
  if(MINGW)
    target_link_libraries(tblastn_gui PRIVATE uuid ole32 oleaut32 imm32 dwmapi shell32)
  endif()
endif()
