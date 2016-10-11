[
@[ for dos in DOS ]@
{
    "property-id" "tag:staff@noreply.openkim.org,2014-05-21:property/phonon-dispersion-dos-cubic-crystal-npt"
    "instance-id" @<dos.iter>@

    "short-name" {
        "source-value"  ["@<crystal_structure>@"]
    }
    "species" {
        "source-value"  ["@<element>@"]
    }
    "a" {
        "source-value"  @<lattice_constant>@
        "source-unit"   "angstrom"
    }
    "basis-atom-coordinates" {
        "source-value"  @<basis_atoms>@
    }
    "space-group" {
        "source-value"  "Fm-3m"
    }
    "temperature" {
        "source-value"  0
        "source-unit"   "K"
    }
    "cauchy-stress" {
        "source-value"  [ 0 0 0 0 0 0 ]
        "source-unit"   "GPa"
    }
    "energy" {
        "source-value"  @<dos.energy|json>@
        "source-unit"   "meV"
    }
    "density-of-states" {
        "source-value"  @<dos.density|json>@
    }
}
@[ endfor ]@

@[ for kpoint in kpoints ]@
{
    "property-id" "tag:staff@noreply.openkim.org,2014-05-21:property/phonon-dispersion-relation-cubic-crystal-npt"
    "instance-id" @<kpoint.iter>@

    "short-name" {
        "source-value"  ["@<crystal_structure>@"]
    }
    "species" {
        "source-value"  ["@<element>@"]
    }
    "a" {
        "source-value"  @<lattice_constant>@
        "source-unit"   "angstrom"
    }
    "basis-atom-coordinates" {
        "source-value"  @<basis_atoms>@
    }
    "space-group" {
        "source-value"  "Fm-3m"
    }
    "temperature" {
        "source-value"  0
        "source-unit"   "K"
    }
    "cauchy-stress" {
        "source-value"  [ 0 0 0 0 0 0 ]
        "source-unit"   "GPa"
    }
    "wave-vector-direction" {
        "source-value" @<wavevector|json>@
        "source-unit"  "1/angstrom"
    }

    "branch-label" {
        "source-value"  ["@<kpoint.branch_label>@"]
    }
    "wave-number" {
        "source-value"  @<kpoint.wavenumber|json>@
        "source-unit"   "1/angstrom"
    }
    "response-frequency" {
        "source-value"  @<kpoint.frequency|json>@
        "source-unit"   "meV"
    }
}
@[ endfor ]@
]