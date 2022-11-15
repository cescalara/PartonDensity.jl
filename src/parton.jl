struct Parton{PT} end
export Parton

const gluon = Parton{:g}()

const d_quark = Parton{:d}()
const d_bar = Parton{:d̄}()
const u_quark = Parton{:u}()
const u_bar = Parton{:ū}()
const s_quark = Parton{:s}()
const s_bar = Parton{:s̄}()
const c_quark = Parton{:c}()
const c_bar = Parton{:c̄}()
const b_quark = Parton{:b}()
const b_bar = Parton{:b̄}()
const t_quark = Parton{:t}()
const t_bar = Parton{:t̄}()

const u_valence = Parton{:uv}()
const d_valence = Parton{:dv}()

export gluon, d_quark, d_bar, u_quark, u_bar, s_quark, s_bar, c_quark, c_bar, b_quark, b_bar, t_quark, t_bar, u_valence, d_valence, sea_quark, u_sea, d_sea


function pdf_map_row end

pdf_map_row(::typeof(gluon)) = 0
pdf_map_row(::typeof(u_valence)) = 1
pdf_map_row(::typeof(d_valence)) = 2
pdf_map_row(::typeof(u_bar)) = 3
pdf_map_row(::typeof(d_bar)) = 4
pdf_map_row(::typeof(s_quark)) = 5
pdf_map_row(::typeof(s_bar)) = 6
pdf_map_row(::typeof(c_quark)) = 7
pdf_map_row(::typeof(c_bar)) = 8
pdf_map_row(::typeof(b_quark)) = 9
pdf_map_row(::typeof(b_bar)) = 10
