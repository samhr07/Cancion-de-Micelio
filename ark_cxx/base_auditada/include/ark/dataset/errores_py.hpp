// errores_py.hpp -- excepciones C++ que llevan el nombre de su clase Python (para el contraste "Clase: mensaje" del porte de
// filtros_de_ventana y deconvolucion_tv). Jefe de Arquitectos, 2026-11-03.
#pragma once
#include <stdexcept>
#include <string>

namespace ark::dataset::ventanas {

struct ErrorPy : std::runtime_error {
    using std::runtime_error::runtime_error;
    virtual const char* clase() const noexcept { return "Exception"; }
};

#define ARK_ERROR_PY(NOMBRE, PADRE)                                              \
    struct NOMBRE : PADRE {                                                      \
        using PADRE::PADRE;                                                      \
        const char* clase() const noexcept override { return #NOMBRE; }          \
    };

ARK_ERROR_PY(TypeError, ErrorPy)
ARK_ERROR_PY(ValueError, ErrorPy)
ARK_ERROR_PY(IndexError, ErrorPy)

}  // namespace ark::dataset::ventanas
