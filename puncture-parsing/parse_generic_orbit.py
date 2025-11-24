from sympy.parsing import mathematica
from sympy import cse, ccode, symbols
from sympy.utilities.iterables import numbered_symbols
from sympy.codegen.rewriting import create_expand_pow_optimization
from collections import deque

with open(f"Kerr/Psi0accKerr.txt") as f:
    psi_str = f.read()
with open(f"Kerr/dtPsi0accKerr.txt") as f:
    dtpsi_str = f.read()
with open(f"Kerr/dxPsi0accKerr.txt") as f:
    dxpsi_str = f.read()
with open(f"Kerr/dyPsi0accKerr.txt") as f:
    dypsi_str = f.read()
with open(f"Kerr/dzPsi0accKerr.txt") as f:
    dzpsi_str = f.read()

# define the symbols that may be needed in the expressions
(
    xp,
    yp,
    zp,
    xpdot,
    xpddot,
    ypdot,
    ypddot,
    zpdot,
    zpddot,
    Dx,
    Dy,
    Dz,
    rp,
    rpdot,
    fx,
    fy,
    ft,
    fxdot,
    fydot,
    ftdot,
    Duft,
    Dufx,
    Dufy,
    Duftdot,
    Dufxdot,
    Dufydot,
    M,
    a,
    l,
) = symbols(
    "xp yp zp xpdot xpddot ypdot ypddot zpdot zpddot Dx Dy Dz rp rpdot fx fy ft fxdot fydot ftdot Duft Dufx Dufy Duftdot Dufxdot Dufydot M a l"
)
psi = mathematica.parse_mathematica(psi_str)
dtpsi = mathematica.parse_mathematica(dtpsi_str)
dxpsi = mathematica.parse_mathematica(dxpsi_str)
dypsi = mathematica.parse_mathematica(dypsi_str)
dzpsi = mathematica.parse_mathematica(dzpsi_str)

subexpr, res = cse([psi, dtpsi, dxpsi, dypsi, dzpsi], optimizations="basic")

cpp_file = """
// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureField.hpp"

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/DynamicBuffer.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Worldtube {{

void acceleration_terms_kerr_0(
    gsl::not_null<Variables<tmpl::list<
        CurvedScalarWave::Tags::Psi, ::Tags::dt<CurvedScalarWave::Tags::Psi>,
        ::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                      Frame::Inertial>>>*>
        result,
    const tnsr::I<DataVector, 3, Frame::Inertial>& centered_coords,
    const tnsr::I<double, 3>& particle_position,
    const tnsr::I<double, 3>& particle_velocity,
    const tnsr::I<double, 3>& particle_acceleration, const double ft,
    const double fx, const double fy, const double fz, const double dt_ft,
    const double dt_fx, const double dt_fy, const double dt_fz,
    const double BH_mass, const std::array<double,3> &BH_spin) {{
  const size_t grid_size = get<0>(centered_coords).size();
  result->initialize(grid_size);
  const double xp = particle_position[0];
  const double yp = particle_position[1];
  const double zp = particle_position[2];
  const double xpdot = particle_velocity[0];
  const double ypdot = particle_velocity[1];
  const double zpdot = particle_velocity[2];

  const double xpddot = particle_acceleration[0];
  const double ypddot = particle_acceleration[1];
  const double zpddot = particle_acceleration[2];

  const double rp = get(magnitude(particle_position));
  const double rpdot = (xp * xpdot + yp * ypdot + zp * zpdot) / rp;

  const auto& Dx = get<0>(centered_coords);
  const auto& Dy = get<1>(centered_coords);
  const auto& Dz = get<2>(centered_coords);

  const double M = BH_mass;
  const double a = BH_spin[2];

    DynamicBuffer<DataVector> temps({number_of_temps}, grid_size);

    {aux_vars}

  get(get<CurvedScalarWave::Tags::Psi>(*result)) = {sol_psi};
  get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(*result)) = {sol_dtpsi};
   get<0>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>, Frame::Inertial>>(
        *result)) = {sol_dxpsi};
   get<1>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>, Frame::Inertial>>(
        *result)) = {sol_dypsi};
   get<2>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>, Frame::Inertial>>(
        *result)) = {sol_dzpsi};
}}
}}

"""
symbol_gen = numbered_symbols()

temp_symbols = {}
for i in range(10000):
    temp_symbols[next(symbol_gen)] = i

temp_dict = {}
root_symbols = [
    Dx,
    Dy,
    Dz,
    xp,
    yp,
    zp,
    xpdot,
    xpddot,
    ypdot,
    ypddot,
    zpdot,
    zpddot,
    rp,
    rpdot,
    ft,
    fx,
    fy,
    fxdot,
    fydot,
    ftdot,
    Duft,
    Dufx,
    Dufy,
    Duftdot,
    Dufxdot,
    Dufydot,
    a,
    M,
]
double_symbols = numbered_symbols("d_")
dv_symbols = numbered_symbols("dv_")
dvs = []
doubles = []
exprs = []
for i, (name, expr) in enumerate(subexpr):
    old_expr = expr
    for sym in expr.free_symbols:
        if sym not in root_symbols:
            expr = expr.subs(sym, exprs[temp_symbols[sym]])
    new_symbol = None
    if Dx in expr.free_symbols or Dy in expr.free_symbols or Dz in expr.free_symbols:
        new_symbol = next(dv_symbols)
        dvs.append((new_symbol, old_expr))
    else:
        new_symbol = next(double_symbols)
        doubles.append((new_symbol, old_expr))
    temp_dict[name] = new_symbol
    exprs.append(expr)

print(f"{len(doubles)} double vars")
print(f"{len(dvs)} datavector vars")

expand_opt = create_expand_pow_optimization(10)
for i, (name, expr) in enumerate(doubles):
    for sym in expr.free_symbols:
        if sym in temp_dict:
            expr = expr.subs(sym, temp_dict[sym])
    doubles[i] = (name, expand_opt(expr))

for i, (name, expr) in enumerate(dvs):
    for sym in expr.free_symbols:
        if sym in temp_dict:
            expr = expr.subs(sym, temp_dict[sym])
    dvs[i] = (name, expand_opt(expr))

for i, expr in enumerate(res):
    for sym in expr.free_symbols:
        if sym in temp_dict:
            expr = expr.subs(sym, temp_dict[sym])
    res[i] = expand_opt(expr)


def check_last_usage(sym, exprs, res):
    for _, expr in exprs:
        if sym in expr.free_symbols:
            return False
    for expr in res:
        if sym in expr.free_symbols:
            return False
    return True


free_storage = deque()
occ_storage = []
total_storage = 0
dv_symbols = [name for name, _ in dvs]
for i, (name, expr) in enumerate(dvs):
    for sym in expr.free_symbols:
        if sym in dv_symbols:
            if check_last_usage(sym, dvs[i + 1 :], res):
                free_storage.append(sym)
    if len(free_storage) > 0:
        occ_storage.append(free_storage.pop())
    else:
        total_storage += 1
        occ_storage.append(name)
print(f"{total_storage} allocations")
free_storage = deque(range(total_storage))
occ_storage = []
current_storage = {}
aux_dvs = ""
for i, (name, expr) in enumerate(dvs):
    for sym in expr.free_symbols:
        if sym in dv_symbols:
            if check_last_usage(sym, dvs[i + 1 :], res):
                free_storage.appendleft(current_storage[sym])
    index = free_storage.popleft()
    current_storage[name] = index
    code = f"DataVector& {name} = temps.at({index});\n{name} = {ccode(expr)};\n"
    aux_dvs += code

aux_doubles = "\n".join(
    f"const double {name} = {ccode(expr)};" for name, expr in doubles
)
aux_vars = aux_doubles + "\n" + aux_dvs

cse_dict = {
    "aux_vars": aux_vars,
    "sol_psi": ccode(res[0]),
    "sol_dtpsi": ccode(res[1]),
    "sol_dxpsi": ccode(res[2]),
    "sol_dypsi": ccode(res[3]),
    "sol_dzpsi": ccode(res[4]),
    "number_of_temps": total_storage,
}

full_file_to_write = cpp_file.format(**cse_dict)
import re

full_file_to_write = re.sub(r"\s(\d)\s", r"\1\.0", full_file_to_write)

with open("Kerr/AccelerationTermsKerrOrder0.cpp", "w") as f:
    f.write(full_file_to_write)
