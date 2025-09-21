# Plan: Convert Documentation Examples to Literate.jl Format

## Objective
Convert the existing markdown documentation examples in `docs/src/examples/` to Literate.jl format to enable executable documentation with automatic code generation and testing.

## Current State
- 5 markdown files in `docs/src/examples/`: `bell.md`, `certify_ground_state_property.md`, `ground_state_energy.md`, `trace_poly.md`, `werner_state.md`
- Current format: Static markdown with embedded Julia code blocks
- Documentation built using Documenter.jl
- No Literate.jl integration currently

## Proposed Approach

### 1. Create Literate.jl Source Files
For each example, create a `.jl` file in a new `examples/literate/` directory with:
- Executable Julia code
- Markdown comments using `# ` for documentation
- Proper Literate.jl structure (code chunks and markdown sections)

### 2. Update Build Configuration
- Add Literate.jl to Project.toml dependencies
- Modify `docs/make.jl` to process Literate.jl files
- Update documentation structure to include generated markdown files

### 3. Conversion Strategy
For each existing markdown file:
- Extract Julia code blocks
- Convert markdown text to Literate.jl comments
- Ensure code is executable and produces expected outputs
- Add proper Literate.jl metadata and structure

## Key Data Structures and Functions

### Directory Structure
```
docs/
├── src/
│   ├── examples/
│   │   ├── literate/     # New: Literate.jl source files
│   │   │   ├── bell.jl
│   │   │   ├── certify_ground_state_property.jl
│   │   │   ├── ground_state_energy.jl
│   │   │   ├── trace_poly.jl
│   │   │   └── werner_state.jl
│   │   └── generated/    # New: Generated markdown files
└── make.jl              # Modified to support Literate.jl
```

### Literate.jl Configuration
- Use default settings with custom preprocessing if needed
- Generate both markdown and notebook outputs
- Maintain existing documentation structure

## Implementation Steps

1. **Setup Literate.jl Integration**
   - Add Literate.jl to Project.toml
   - Create directory structure
   - Update make.jl with Literate.jl processing

2. **Convert First Example (bell.md)**
   - Create bell.jl with Literate.jl format
   - Test code execution
   - Generate markdown and verify output

3. **Convert Remaining Examples**
   - Convert certify_ground_state_property.jl
   - Convert ground_state_energy.jl
   - Convert trace_poly.jl
   - Convert werner_state.jl

4. **Update Documentation Build**
   - Modify make.jl to process all Literate.jl files
   - Test full documentation build
   - Verify all examples execute correctly

5. **Cleanup and Validation**
   - Remove old static markdown files
   - Update documentation navigation
   - Test final documentation build

## Expected Outcomes
- Executable documentation examples
- Automatic code output verification
- Improved maintainability
- Better testing of documentation code