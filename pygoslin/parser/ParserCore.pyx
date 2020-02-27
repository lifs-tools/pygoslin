
def parser_core(text_to_parse, dict nt, dict t, list lft, list rgt):
    
    cdef:
        list DP, ld, Ks, DL, DR
        dict lld
        set s
        long rule_index, mask, new_key, old_key, i_shift, index_pair_1, index_pair_2, key
        int n = len(text_to_parse), i, j, k, jpok, im1, shift
    
    Ks = [set([0]) for i in range(n)]
    DP = [[{} for j in range(n - i)] for i in range(n)]
    DL = [[set() for j in range(n - i)] for i in range(n)]
    DR = [[set() for j in range(n - i)] for i in range(n)]
    
    shift, mask = 32, (1 << 32) - 1    
    for i in range(n):
        c = text_to_parse[i]
        if c not in t: return None
        
        for rule_index in t[c]:
            new_key = rule_index >> 32
            old_key = rule_index & mask
            dp_node = [c, None, old_key, None]
            DP[i][0][new_key] = dp_node
            DL[i][0] |= lft[new_key]
            DR[i][0] |= rgt[new_key]
    
    
    for i in range (1, n):
        im1 = i - 1
        for j in range(n - i):
            
            DPj, jp1 = DP[j], j + 1
            DPji = DPj[i]
            Ksj = Ks[j]
            
            for k in Ksj:
                jpok = jp1 + k
                if im1 - k not in Ks[jpok]: continue
            
                D1, D2 = DPj[k], DP[jpok][im1 - k]
                
                for key in DR[j][k] & DL[jpok][im1 - k]:
                    index_pair_1 = key >> shift
                    index_pair_2 = key & mask
                    parse_content = [index_pair_1, D1[index_pair_1], index_pair_2, D2[index_pair_2]]
                    for rule_index in nt[key]:
                        DPji[rule_index] = parse_content
                        DL[j][i] |= lft[rule_index]
                        DR[j][i] |= rgt[rule_index]
                                    
            if DPji: Ksj.add(i)
            
            
    return DP