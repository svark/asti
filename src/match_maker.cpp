#include <list>
#include <limits>
// given a set of source vectors(svecs) and target vectors (tvecs) that
// represent discrete segments between corners in source and in target:
// with number of source vectors being more than in target, find out
// indices of matched pairs (i,j) such that partial sums of the source
// and target vectors upto i and j resp, can be matched with minimum
// distortion. Distortion is a measure of how disoriented the segments
// are
#if 0
template<class VectorContT>
struct vector_match_maker {
    vector_match_maker( const VectorContT & _svecs,
                        const VectorContT & _tvecs )
        :src_vec_diffs(_svecs ),
         tgt_vec_diffs(_tvecs )
    {
        svec_cum_sums = src_vec_diffs;
        tvec_cum_sums = tgt_vec_diffs;
        for(int i = 1; i < src_vec_diffs.size(); ++i)
        {
            svec_cum_sums[i] += svec_cum_sums[i-1];
        }
        for(int j = 1; j < tgt_vec_diffs.size(); ++j)
        {
            tvec_cum_sums[j] += tvec_cum_sums[j-1];
        }
        setup_length_distortion();
    }

    struct matched_pair {
        matched_pair() {
            _i = 0; _j = 0;
            _distortion = -1;
            _next =0; _prev = 0;
        }
        matched_pair( int i , int j ):
            _i(i), _j(j)
        {
            _distortion = -1;
            _next = 0; _prev = 0;
        }
        int _i;
        int _j;
        double _distortion;
    };

    typedef std::list<matched_pair> matched_pairs_t;
    typedef std::list<matched_pair>::iterator matched_pairs_it;

    double distortion_so_far( matched_pairs_it mpr )
    {
        double distortion  = 0;
        for(;mpr != matched_pairs.begin(); --mpr)
        {
            distortion += std::prev(mpr)->_distortion;
        }
        return distortion;
    }

    bool
    update_min_distortion( matched_pairs_it mpr,
                           double &minDistortion,
                           size_t *bestI
        )
    {
        if( mpr->_j != tgt_vec_diffs.size() - 1 )
            return false ;

        if( mpr->_distortion < 0.0 )
            return false ; // invalid path

        double distortion = mpr->_distortion;
        distortion += distortion_so_far(mpr);

        if(minDistortion > distortion )
        {
            minDistortion = distortion;
            size_t *bi = bestI;
            for(auto mp : matched_pairs)
            {
                *bi++=  mp.i;
            };
        }
        return true ;
    }

    // try a different path if we hit a wall.
    matched_pair_t
    try_another_prev_pair(int max_ij_diff, matched_pairs_it mpr)
    {
        while( mpr != matched_pairs.end()
               && mpr->_i - mpr->_j > max_ij_diff )
        {
            // done with previous track, now go back and try a
            // different route
            if(mpr != matched_pairs.begin()) {
                --mpr;
                matched_pairs.erase(std::next(mpr), matched_pairs.end());
                ++mpr->_i;
                mpr->_distortion = -1;
            }else {
                matched_pairs.erase(mpr, matched_pairs.end());
                mpr = matched_pairs.end();
            }
        }
        return mpr;
    }

    // compute  vector sums starting from last matched source corner.
    // up until  corner at mpr->_i-1
    vec_t
    setup_base_vec(matched_pair_it mpr)
    {
        assert(mpr!= matched_pairs.end());
        vec_t base_src_vec(0.0);
        if( mpr->_i == 0 )
            return base_src_vec;

        int start_i =
            ((mpr !=
              matched_pairs.begin())? (std::prev(mpr->_i) + 1): 0);

        if(mpr->_i == start_i) {
            return base_src_vec;
        }

        base_src_vec = svec_cum_sums[mpr->_i-1];

        if(start_i > 0){
            base_src_vec = base_src_vec - svec_cum_sums[start_i-1] ;
        }

        return base_src_vec;
    }

    // find current distortion at the i,j corner and set up a new
    // context for matching the next corner at target at index j+1

    matched_pairs_it
    update_distortion(int max_ij_diff, matched_pairs_it mpr,
                      double minDistortion)
    {
        vec_t total_tvec (tvec_cum_sums.back());
        auto total_tvec_len =  len(total_tvec);
        int j = mpr->_j;
        vec_t tvec (tgt_vec_diffs[j] );
        // ensure last corner in target  is matched with last in source
        int &i = mpr->_i; // NOTE the reference
        if( j == tgt_vec_diffs.size() - 1  &&
            i < src_vec_diffs.size() - 1 )
            i = src_vec_diffs.size() - 1;

        vec_t total_svec(setup_base_vec ( mpr ));
        double distSoFar  = distortion_so_far( mpr );
        double tlen       = len(tvec);
        double denom      = tlen * total_tvec_len ;
        double cos_theta2 = dot(tvec, total_tvec)/denom;
        double sintheta2  = sqrt(1 - cos_theta2 * cos_theta2);
        for( ; i - j <= max_ij_diff; ++i )
        {
            const vec_t &svec = src_vec_diffs[i];
            total_svec = total_svec + svec;
            double slen = len(total_svec);

            if( slen > len_distortion*tlen ||
                tlen > len_distortion*slen ) {
                // dont try too hard..
                continue;
            }
            mpr->_distortion = -1;
            double cos_theta1 = dot(total_svec, total_tvec)
                / (slen * total_tvec_len);

            double sin_theta1 = sqrt(1 - cos_theta1 * cos_theta1);
            double cos_theta2_theta1 = cos_theta2 * cos_theta1
                + sin_theta1 * sintheta2;

            double cur_distortion = (1 - cos_theta2_theta1 );

            if(cur_distortion + distSoFar > minDistortion) {
                continue;
            }

            mpr->_distortion = cur_distortion;
        }
        return mpr;
    }

    // length distortion may occur when source and target distances
    // between corners are not same. To match corners we assume that
    // length distortion is not much. By default we assume that upto
    // distances are within 1.5 times of each other. Using too large
    // values slows down the search for matching pairs of corners and
    // too low would mean losing out on accuracy in matching. But to
    // guess a value we just look at the distance to first and last but
    // one corner in target and min and max lengths in the source
    void setup_length_distortion()
    {
        if(src_vec_diffs.size() == 0 || tgt_vec_diffs.size() == 0)
            return;
        double max_src_dist = 0;
        double min_src_dist = len(src_vec_diffs[0]);
        ptrdiff_t max_ij_diff =
            src_vec_diffs.size() -
            tgt_vec_diffs.size();

        len_distortion = 1.3;

        for( int i = 1 ; i <= max_ij_diff; ++i )
        {
            double len = len(svec_cum_sums[i]);
            max_src_dist = std::max(max_src_dist, len);
            min_src_dist = std::min(min_src_dist, len);
        }

        double ftgt_len =  len(tgt_vec_diffs[0]);
        if(len_distortion < ftgt_len / max_src_dist )
            len_distortion = ftgt_len / max_src_dist;

        if( len_distortion < min_src_dist/ ftgt_len )
            len_distortion = min_src_dist/ ftgt_len ;

        max_src_dist = 0;
        size_t l = src_vec_diffs.size() -1;
        min_src_dist = len(src_vec_diffs[l]);

        for( int i = 1 ; i <= max_ij_diff; ++i )
        {
            double len_ = len(svec_cum_sums[l-i] -
                              svec_cum_sums[l] );
            max_src_dist = std::max(max_src_dist, len_ );
            min_src_dist = std::min(min_src_dist, len_ );
        }

        size_t  m = tgt_vec_diffs.size()-1;
        auto tgt_vec_last_len = len(tgt_vec_diffs[m]);
        if(len_distortion < tgt_vec_last_len / max_src_dist  )
            len_distortion = tgt_vec_last_len / max_src_dist;

        if( len_distortion < min_src_dist/ tgt_vec_last_len )
            len_distortion = min_src_dist/ tgt_vec_last_len;
    }

    // Having matched a pair (i,j) start searching for next matching
    // source corner for j+1 corner in target, starting from i+1

    matched_pair_it
    try_next_pair(int max_ij_diff, matched_pair_it mpr)
    {
        if(mpr != matched_pairs.end() &&
           mpr->_i - mpr->_j <= max_ij_diff) {
            // havent hit a wall... so try next pair
            int i = mpr->_i;
            int j = mpr->_j;
            if( j < tgt_vec_diffs.size()-1)
            {
                matched_pairs.insert(std::next(mpr),
                                     new matched_pair(i+1, j+1));
                ++mpr;
            }
            else
            {
                // this is to force the matchmaker to try another path...
                ++mpr->_i;
                mpr->_distortion = -1;
            }
        }
        return mpr;
    }

    // find best (i,j) pairs (i--> source vector index, j-> target vec
    // index) that minimizes overall distortion in orientation of
    // matched vector pairs. Distortion is taken as sum of = 1 - cosine
    // of angle between the two source and target vectors between
    // matching corners.
    double match_vectors(size_t *bestI,
                         double minTotalDistortion  )
    {
        matched_pairs =  matched_pairs_t();
        matched_pairs.push_back(matched_pair());
        if(src_vec_diffs.size() == 0 ||
           tgt_vec_diffs.size() == 0)
            return std::numeric_limits < double >::infinity;

        int max_ij_diff = src_vec_diffs.size() - tgt_vec_diffs.size();

        // assuming i matches j,
        double minDistortion = minTotalDistortion;

        auto mpr =  matched_pair.begin();
        // big for loop to avoid recursively calling
        // match_vectors..
        while(mpr != matched_pair.end() ) {

            // update current distortion and match next corner at target
            mpr = update_distortion(max_ij_diff, mpr,
                                   minDistortion);

            // have we reached a full match with all corners
            // matching. If yes update minDistortion..
            update_min_distortion(mpr, minDistortion, bestI);

            mpr = try_next_pair(max_ij_diff, mpr);

            // have we exhausted all possibilties..
            // then go back and try a new a matched pair..
            mpr = try_another_prev_pair(max_ij_diff, mpr);
        }
        len_distortion += 0.5;
        return minDistortion;
    }
private:
    const VectorContT &src_vec_diffs;
    const VectorContT &tgt_vec_diffs;
    VectorContT svec_cum_sums;
    VectorContT tvec_cum_sums;
    double len_distortion;
    matched_pairs_t matched_pairs;
};

#endif
