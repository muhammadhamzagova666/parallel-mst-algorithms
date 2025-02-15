# Parallel MST Algorithms  
*Efficient Parallel Implementations of Boruvka's & Kruskal's Algorithms*

## Overview
This project implements and compares parallel algorithms for computing the Minimum Spanning Tree (MST) of weighted graphs. Leveraging OpenMP and MPI, the repository features both sequential and parallel versions of Boruvka's and Kruskal's algorithms. Designed for developers, researchers, and academics, this project highlights scalability, performance metrics, and design trade-offs of parallel graph processing.

**Key Features & Functionalities:**
- **Parallel Processing:** Uses OpenMP and MPI to optimize graph processing tasks.
- **Multiple Algorithm Implementations:** Provides sequential and parallel versions of Boruvka's and Kruskal's MST algorithms.
- **Performance Evaluation:** Analyzes execution, computation, and communication times across different dataset sizes and processor counts.
- **Comprehensive Documentation:** Detailed documentation with performance insights, algorithmic background, and code architecture.

**Target Audience:**  
Developers, researchers, and students interested in parallel algorithms, graph theory, and performance optimization.

**Unique Selling Points:**  
- In-depth comparison of parallel versus sequential implementations.
- Real-world performance metrics and empirical evaluations.
- Clean, modular source code for easy adaptation and enhancement.

## Technology Stack
- **C Language** (Standard C99)
- **OpenMP** (Parallelization)
- **MPI** (Message Passing Interface for distributed computing)
- **LaTeX** (Documentation, report generation)
- **GNU Tools** (Compilation and debugging)
- **Others:** Standard libraries for memory management and file I/O.

## Installation & Setup
### Prerequisites
- GNU Compiler Collection (GCC)
- OpenMP compatible compiler (e.g., `gcc`)
- MPI library (e.g., `mpich` or `openmpi`) for MPI implementations
- LaTeX distribution for report generation (e.g., TeX Live)

### Steps to Install
1. **Clone the Repository:**
   ```sh
   git clone https://github.com/muhammadhamzagova666/parallel-mst-algorithms.git
   cd parallel-mst-algorithms
   ```

2. **Build the Project:**
   For example, to compile the OpenMP version of Kruskal's algorithm:
   ```sh
   gcc -fopenmp Source/kruskal_omp.c Source/Graph.c -o kruskal_omp
   ```
   To compile the MPI version of Kruskal's algorithm:
   ```sh
   mpicc Source/kruskal_mpi.c Source/Graph.c -o kruskal_mpi
   ```

3. **Generate the Report:**
   Navigate to the Report directory and compile the LaTeX file:
   ```sh
   cd Report
   pdflatex PDC.tex
   ```

## Usage Guide
- **Running Algorithms:**
  - For the **OpenMP** version:
    ```sh
    ./kruskal_omp
    ```
    Follow on-screen prompts to input graph node count.

  - For the **MPI** version:
    ```sh
    mpirun -np 4 ./kruskal_mpi
    ```
    The number of processes (`-np`) can be varied as needed.

- **Example Output:**
  - Total MST cost and execution time are printed.
  - Graphical outputs and performance metrics are provided in the generated report.

## Project Structure
```
parallel-mst-algorithms/
├── Parallel Versions of Boruvka’s and Kruskal's Algorithms Proposal.docx
├── Parallel Versions of Boruvka’s and Kruskal's Algorithms Report.pdf
├── Parallel Versions of Boruvka’s and Kruskal's AlgorithmsGraphs.xlsx
├── Report/
│   ├── boruvka_algo.PNG
│   ├── boruvka_algo2.PNG
│   ├── boruvka_graph.jpeg
│   ├── FAST.png
│   ├── kruskal_algo.PNG
│   ├── kruskal_graph.jpeg
│   ├── mst_algo.png
│   ├── NU-logo.jpg
│   ├── PDC.tex
│   └── sample.png
└── Source/
    ├── boruvka_openmp.c
    ├── boruvka_s.c
    ├── boruvka.c
    ├── Graph.c
    ├── kruskal_mpi.c
    ├── kruskal_omp.c
    └── kruskal_s.c
```
**Key Directories and Files:**
- **Source:** Contains source code for both sequential and parallel algorithm implementations.
- **Report:** Contains project documentation, images, and the LaTeX report file.
- **Documentation Files:** Detailed proposal and report documents outlining algorithmic approach and performance analysis.

## Configuration & Environment Variables
- **.env / config.json:** (If applicable in future versions for dynamic graph data, etc.)
- **Example Configuration:**
  ```json
  {
    "num_nodes": 1000,
    "density": 0.5,
    "max_threads": 8
  }
  ```

## Deployment Guide
- **Local Deployment:**
  Run the compiled executables locally on your workstation.
- **Containerized Deployment:**
  A Dockerfile can be added to containerize the project for easier deployment.
- **CI/CD Integration:**
  Integration with GitHub Actions for automated builds and testing is recommended.

## Testing & Debugging
- **Running Tests:**
  - Unit tests can be integrated using frameworks like CUnit.
  - Example command:
    ```sh
    make test
    ```
- **Debugging Tips:**
  - Use `gdb` for stepping through code.
  - Print statements and logging are included to track performance metrics.

## Performance Optimization
- Use profiling tools (e.g., `gprof`) to identify bottlenecks.
- Experiment with different process counts and thread numbers.
- Optimize memory allocations and reduce redundant operations.

## Security Best Practices
- Validate all user inputs.
- Ensure proper memory management to avoid leaks and buffer overflows.
- Regularly update dependencies.

## Contributing Guidelines
- **How to Contribute:**
  - Fork the repository.
  - Create a feature branch (`git checkout -b feature/your-feature`).
  - Commit your changes following clear commit messages.
  - Open a pull request for review.
- **Issue Reporting:**
  - Please open issues for bugs, feature requests, or improvements.
- **Code of Conduct:**
  - Maintain professionalism and respect in all interactions.

## Documentation
- Detailed documentation can be found in the `docs/` directory (to be added).
- API references and additional design details are provided in the report.

## Roadmap
- Integration of automated tests.
- Docker container support.
- Enhancement of performance profiling and scalability tests.
- Extended documentation and user guides.

## FAQ
**Q:** What are the prerequisites for running the project?  
**A:** You need a C compiler with OpenMP support, an MPI library, and a LaTeX distribution to render the report.

**Q:** How can I contribute to this project?  
**A:** Please see the Contributing Guidelines section above.

## Acknowledgments & Credits
- Special thanks to all contributors and the open-source community.
- Libraries and tools such as OpenMP, MPI, and LaTeX are gratefully acknowledged.

## Contact Information
For any inquiries or support, please reach out via GitHub or contact:
- GitHub: [muhammadhamzagova666](https://github.com/muhammadhamzagova666)

---

*Happy coding and exploring parallel graph algorithms!*
