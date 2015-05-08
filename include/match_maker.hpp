namespace geom
{

template <class VectorContT>
extern std::vector<std::pair<size_t, size_t > >
match_vectors(const VectorContT & vs1, const VectorContT & vs2 );

template <class PointContT>
extern std::vector<std::pair<size_t, size_t > >
match_points(const PointContT & ps1, const PointContT & ps2 )
{
    std::vector<RAWTYPE(ps1[0] -  ps1[1]) > vs1(ps1.size());
    std::vector<RAWTYPE(ps2[0] -  ps2[1]) > vs2(ps2.size());
    std::adjacent_difference(ps1.begin(), ps1.end(), vs1.begin());
    std::adjacent_difference(ps2.begin(), ps2.end(), vs2.begin());
    return match_vectors(util::mk_range(vs1.begin() + 1, vs1.end()),
                         util::mk_range(vs2.begin() + 1, vs2.end()));
}

}
