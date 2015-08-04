namespace lu
{
template <typename MatrixT>
void
banded_lu_decompose(MatrixT & _A)
{
    size_t n = _A.rows();
    for(size_t k = 0; k < n - 1; ++k) {
        _A(k + 1, k)     = _A(k + 1, k) / _A(k, k);
        _A(k + 1, k + 1) -=  _A(k + 1, k) * _A(k, k + 1);
    }
}

template <typename MatrixT, typename VectorT>
void
banded_lu_solve(const MatrixT & _A, VectorT & b)
{
	size_t n = _A.rows();
    for(size_t j = 0;j < n; ++j)
    {
        b[j + 1] -= _A(j + 1, j) * b[j];
    }
    for(size_t j = n - 1;j >= 0; --j)
    {
        b[j] = b[j] / _A(j, j) ;
        b[j - 1] -= _A(j - 1, j) * b[j];
    }
    return ;
}
}
