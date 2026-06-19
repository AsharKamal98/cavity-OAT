# Parser

This file gives instructions for creating project data classes in parser files.
These classes are used for data validation, clear function outputs, and easy
attribute-style access to computed data.

The parsing file can be found in `docs/instructions/parser.py`. If one does not
already exist, create one.

## Rules

- Use Pydantic classes for new parser output containers.
- Use explicit field types for all saved data.
- Use the nullable union convention for optional fields.
- Use `= None` when the default value is `None`.
- Prefer named fields over generic dictionaries or lists when the output shape
  is part of the public API.
- Use consistent class names that make related classes easy to recognize.

Follow this style for future parser output classes:

```python
class SMomentSeries(BaseModel):
    phase_index: Array
    Sx: Array
    Sy: Array


class JMomentSeries(BaseModel):
    phase_index: Array
    Jx: Array
    Jy: Array


class MomentSeries(BaseModel):
    t: Array
    J: JMomentSeries | None = None
    S: SMomentSeries | None = None
```
