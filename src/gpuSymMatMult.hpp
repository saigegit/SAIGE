/*! @file gpuSymMatMult.hpp
 *  @author Youngdae Kim
 *  @brief Perform distributed matrix operations
 */

#ifndef GPU_SYM_MAT_MULT_HPP
#define GPU_SYM_MAT_MULT_HPP

/*! Computes G*x where G is a GRM formed by G=AA^T by A*(A^T*x) using cublas on GPUs.
 *  The class gpuSymMatMult stores a submatrix of A and provides a routine sym_sgemv for G*x.
 *  The routine sym_sgemv_range is for submatrix multiplications.
 */
class gpuSymMatMult
{
public:
    bool m_isLoaded;           // true if data is loaded in GPU memory
    size_t m_rows;             // number of rows of a submatrix of A
    size_t m_cols;             // number of columns of a submatrix of A
    size_t m_global_cols;      // number of columns of A
    size_t m_global_col_start; // start global column index of a submatrix of A
    void *m_handle;            // handle for cublas
    float *m_A;                // pointer to data of this matrix
    float *m_x;                // buffer for sgemv
    float *m_y;                // buffer for sgemv
    float *m_z;                // buffer for sgemv

    gpuSymMatMult();
    ~gpuSymMatMult();

    bool is_loaded() { return m_isLoaded; }
    int set_matrix(int rank, size_t g_col_start, size_t n_rows, size_t n_cols, const float *A);
    int sym_sgemv(int rank, size_t n_elem, const float *x, float *ret);
    int sym_sgemv_range(int rank, size_t g_start_col, size_t g_end_col, size_t n_elem, const float *x, float *ret);
};
#endif // GPU_SYM_MAT_MULT_HPP
