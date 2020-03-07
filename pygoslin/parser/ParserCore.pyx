# distutils: language=c++

cimport cython

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)
def parser_core(unicode text_to_parse, dict nt, dict t, list lft, list rgt):
    
    cdef:
        list DP, ld, Ks, DL, DR
        dict lld
        long rule_index, mask, index_pair_1, index_pair_2, key
        int n = len(text_to_parse), i, j, k, jpok, im1, shift
        unicode c
    
    Ks = [set([0]) for i in range(n)]
    DP = [[{} for j in range(n - i)] for i in range(n)]
    DL = [[set() for j in range(n - i)] for i in range(n)]
    DR = [[set() for j in range(n - i)] for i in range(n)]
    
    
    shift, mask = 32, (1 << 32) - 1
    for i, c in enumerate(text_to_parse):
        if c not in t: return None
        
        for rule_index in t[c]:
            dp_node = [c, None, rule_index, None]
            DP[0][i][rule_index] = dp_node
            DL[0][i] |= lft[rule_index]
            DR[0][i] |= rgt[rule_index]
    
    
    for i in range (1, n):
        im1, DPi = i - 1, DP[i]
        for j in range(n - i):
            jp1, DPij, Ksj, DLij, DRij = j + 1, DPi[j], Ks[j], DL[i][j], DR[i][j]
            
            for k in Ksj:
                jpok, im1mk = jp1 + k, im1 - k
                if im1mk not in Ks[jpok]: continue
                
                for key in DR[k][j] & DL[im1mk][jpok]:
                    index_pair_1 = key >> shift
                    index_pair_2 = key & mask
                    parse_content = [index_pair_1, DP[k][j][index_pair_1], index_pair_2, DP[im1mk][jpok][index_pair_2]]
                    for rule_index in nt[key]:
                        DPij[rule_index] = parse_content
                        DLij |= lft[rule_index]
                        DRij |= rgt[rule_index]
                                
            if DPij: Ksj.add(i)

    return DP
