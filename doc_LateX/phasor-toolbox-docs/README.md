# Phasor Toolbox Documentation

This project contains the documentation for the Phasor Toolbox, which implements the `PhasorSS` and `dSpaceDataExplorer` classes. The documentation is structured into chapters, each focusing on different aspects of the toolbox.

## Project Structure

- **src/chapters/introduction.tex**: Introduction to the toolbox, outlining its purpose and scope.
- **src/chapters/phasorSS.tex**: Detailed documentation of the `PhasorSS` class, including properties, methods, and usage examples.
- **src/chapters/dSpaceDataExplorer.tex**: Documentation for the `dSpaceDataExplorer` class, covering its properties, methods, and usage.
- **src/main.tex**: The main LaTeX file that compiles the entire documentation.
- **bibliography.bib**: Bibliography entries formatted in BibTeX for citations within the documentation.

## Building the Documentation

To build the documentation, ensure you have a LaTeX distribution installed (e.g., TeX Live, MiKTeX). Then, navigate to the `src` directory and run the following command:

```bash
pdflatex main.tex
```

This will generate a PDF document containing the complete documentation for the Phasor Toolbox.

## Dependencies

Make sure to have the following packages installed in your LaTeX distribution:

- `amsmath`
- `graphicx`
- `hyperref`
- `biblatex`

These packages are necessary for proper formatting and referencing within the documentation.

## Contributing

If you would like to contribute to the documentation or the toolbox itself, please fork the repository and submit a pull request with your changes.