# Parser

This file gives instructions for creating project data classes in parser files.
These classes are used for data validation, clear function outputs, and easy
attribute-style access to computed data.

Parser classes should live in the repository's parser module or package. In
this repo, parser classes live under `parser/` and are split by domain, for
example `parser/common.py`, `parser/moments.py`, `parser/j_moments.py`,
`parser/mfe.py`, and `parser/mfe_residuals.py`.

## Rules

- Use Pydantic classes for new parser output containers.
- Use explicit field types for all saved data.
- Use the nullable union convention for optional fields.
- Use `= None` when the default value is `None`.
- Use modern built-in generic types such as `list[T]`, `tuple[T, ...]`,
  `dict[K, V]`, and `T | None`.
- Prefer named fields over generic dictionaries or lists when the output shape
  is part of the public API.
- Use consistent class names that make related classes easy to recognize.

Follow this style for future parser output classes:

```python
class JMomentSeries(BaseModel):
    t: Array
    phase_index: Array
    x: Array
    y: Array
    z: Array
    x_groups: tuple[Array, ...] | None = None
    y_groups: tuple[Array, ...] | None = None
    z_groups: tuple[Array, ...] | None = None


class MomentSeries(BaseModel):
    t: Array
    parameters: MomentParameters | None = None
    J: JMomentSeries | None = None
    S: Any | None = None
```