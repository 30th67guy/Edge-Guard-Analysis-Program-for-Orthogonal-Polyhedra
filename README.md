# Edge Guard Analysis Program for Orthogonal Polyhedra

This repository contains the implementation of the Edge Guarding algorithm for orthogonal polyhedra, as presented in the thesis "Statistical Analysis on the Art Gallery Problem for Orthogonal Polyhedra".

In src directory, there are 2 files. The file named as "main_with_viewer.cpp" is the older file.
It is useful to test changing some parameters because it allows us to see a generated polyhedra through GLFW.

## Requirements

- Windows 10 or 11 (x64)

- Visual Studio Code
- Git
- Microsoft C++ Build Tools (MSVC v143)
- CMake
- vcpkg (for library management)

And the following libraries are required.
They are managed via vcpkg.

- GLFW3
- CGAL (Computational Geometry Algorithms Library)

## Build Instructions

To reproduce the environment and build the project, please follow these steps:

### 1. Prerequisites Installation

Visual Studio Code: Install and choose "Add to PATH".

Microsoft C++ Build Tools: Select "Desktop development with C++" and then ensure "MSVC v143", "Windows SDK", and "C++ CMake" are selected.

Git: Install and choose "Git from the command line and also from 3rd-party software".

### 2. Setup vcpkg and Libraries

Run the following commands in your command prompt (e.g., in `C:/Libraries`):

git clone [https://github.com/microsoft/vcpkg.git](https://github.com/microsoft/vcpkg.git)
cd vcpkg
bootstrap-vcpkg.bat
vcpkg integrate install

vcpkg install glfw3 --triplet x64-windows-release
vcpkg install cgal --triplet x64-windows-release

### 3. Project Setup

Clone this repository or download the source code.

Create a build folder inside the project directory.

Then, copy 3 files, "glfw3.dll", "gmp-10.dll" and "mpfr-6.dll", from C:\Libraries\vcpkg\installed\x64-windows-release\bin\ (vcpkg installation folder) to the build folder.

### 4. Build with VS Code

Open the project folder in VS Code.

Install extensions: C/C++ (ms-vscode.cpptools) and CMake Tools (ms-vscode.cmake-tools).

Press Ctrl+Shift+P and select CMake: Select a Kit -> Choose Visual Studio Build Tools 2022 Release - amd64.

Press F7 (or CMake: Build) to compile.

Press F5 (or CMake: Run Without Debugging) to run the program.

### License

This software is released under the GNU General Public License v3.0 (GPLv3). See the LICENSE file for details.
