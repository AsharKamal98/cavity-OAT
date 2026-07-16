#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Phase Protocol]
]

= Purpose

This file specifies the shared phase-protocol hierarchy across solvers
and plotting. The implementation belongs in `parser/common.py` and
`common/utils/phases.py`.

= Phase Hierarchy

`Phase` is one constant integration segment. A `FamilyPhase` is one physical
phase shown in plots and contains one or more `Phase` objects. A
`PhaseProtocol` contains the ordered `FamilyPhase` objects:

```python
Phase(duration, omega, delta, label="")

FamilyPhase(
    duration,
    omega,
    delta,
    ramp_duration=0,
    num_ramp_segments=0,
    label="",
    segments: tuple[Phase, ...] = (),
)

PhaseProtocol(
    family_phases: tuple[FamilyPhase, ...],
)
```

`FamilyPhase.omega` and `FamilyPhase.delta` are its target control values, and
its duration is the sum of its integration `Phase` durations. The protocol
provides the flattened `integration_phases` and the matching
`integration_to_family_index` mapping, together with `total_duration`. When no
ramp is used, `FamilyPhase.segments` contains exactly one target-value `Phase`.

= Method in Pseudo-code

```python
class PhaseProtocol(BaseModel):
    family_phases: tuple[FamilyPhase, ...]

    @staticmethod
    def _build_segments(
        start_values, target_values, duration,
        ramp_duration, num_ramp_segments, label,
    ):
        if num_ramp_segments == 0:
            return (Phase(duration, *target_values, label),)
        ramp_segments = tuple(
            Phase(...)
            for segment_values in interpolated_values  # Includes target values.
        )
        return ramp_segments + optional_target_hold_segment

    @model_validator(mode="after")
    def build_family_segments(self) -> "PhaseProtocol":
        start_values = (0, 0)
        completed_family_phases = []
        for family_phase in self.family_phases:
            target_values = (family_phase.omega, family_phase.delta)
            segments = self._build_segments(
                start_values,
                target_values,
                family_phase.duration,
                family_phase.ramp_duration,
                family_phase.num_ramp_segments,
                family_phase.label,
            )
            completed_family_phases.append(
                family_phase.model_copy(update={"segments": segments})
            )
            start_values = target_values
        self.family_phases = tuple(completed_family_phases)
        return self

def default_three_phase_protocol(
    durations, ..., ramp_durations=(0, 0, 0), ramp_segment_counts=(0, 0, 0)
) -> PhaseProtocol:
    return PhaseProtocol(family_phases=(
        FamilyPhase(durations[0], Omega0, 0, ramp_durations[0], ramp_segment_counts[0]),
        FamilyPhase(durations[1], Omega0, delta0, ramp_durations[1], ramp_segment_counts[1]),
        FamilyPhase(durations[2], 0, 0, ramp_durations[2], ramp_segment_counts[2]),
    ))
```

`PhaseProtocol._build_segments(...)` owns the complete segment-construction
logic. The model validator calls it automatically for every `FamilyPhase`. The
final ramp segment always contains the target values. If `ramp_duration` is
shorter than the family-phase duration, one target-value hold segment fills the
remaining time. With zero ramp duration and zero ramp segments, one target-value
segment spans the full phase. The first `FamilyPhase` starts from $(0, 0)$;
each later `FamilyPhase` starts from the previous `FamilyPhase` target.

= Use by Backend

- `SimulationMetadata` and solver parameter containers receive and store the
  complete protocol as the required `phase_protocol` field.
- MCWF, MFE, and QuTiP evolve through the flattened sequence formed from each
  `FamilyPhase.segments` tuple, exposed as
  `phase_protocol.integration_phases`.
- MCWF precomputes and evolves one integration `Phase` at a time. MFE runs one
  adaptive solve and selects the active integration `Phase` from time. QuTiP
  includes all integration `Phase` boundaries in its solver time grid.
- Plotting and `FamilyPhase`-level summaries consume the `FamilyPhase` objects
  directly through `phase_protocol.family_phases`.
- `MomentSeries` uses `phase_protocol.total_duration` when constructing its
  saved-time grid.
- `phase_boundary_times(...)` returns cumulative boundaries for the supplied
  `FamilyPhase` or integration `Phase` sequence, while
  `phase_values_at_time(...)` evaluates integration `Phase` objects.
- `integration_phase_indices_at_times(...)` maps saved times to their
  `integration_phase_index` values.

= Invariants

- Each integration `Phase` is piecewise constant.
- `ramp_duration` must lie between zero and the family-phase duration.
- `ramp_duration` and `num_ramp_segments` must either both be zero or both be
  positive.
- All backends use the same integration `Phase` sequence for one run.
- Plots remain organized by `FamilyPhase` objects, normally the three physical
  phases.
- A `FamilyPhase` without a ramp has one integration `Phase` and retains the
  current fast precomputation path.
