

def parser_core(text_to_parse, dict nt, dict t, list lft, list rgt, dict substitution):
    
    cdef:
        list DP, ld, Ks, DL, DR
        dict lld
        long rule_index, mask, index_pair_1, index_pair_2, key, s
        int n = len(text_to_parse), i, j, k, jpok, im1, shift
    
    Ks = [set([0]) for i in range(n)]
    DP = [[{} for j in range(n - i)] for i in range(n)]
    DL = [[set() for j in range(n - i)] for i in range(n)]
    DR = [[set() for j in range(n - i)] for i in range(n)]
    
    shift, mask = 32, (1 << 32) - 1    
    for i in range(n):
        c = text_to_parse[i]
        if c not in t: return None
        
        rule_index = t[c]
        dp_node = [c, None, rule_index, None]
        DP[i][0][rule_index] = dp_node
        DL[i][0] |= lft[rule_index]
        DR[i][0] |= rgt[rule_index]
    
        if rule_index in substitution:
            for s in substitution[rule_index]:
                DP[i][0][s] = dp_node
                DL[i][0] |= lft[s]
                DR[i][0] |= rgt[s]
    
    
    for i in range (1, n):
        im1 = i - 1
        for j in range(n - i):
            
            DPj, jp1 = DP[j], j + 1
            DPji, Ksj = DPj[i], Ks[j]
            
            for k in Ksj:
                jpok = jp1 + k
                if im1 - k not in Ks[jpok]: continue
                
                for key in DR[j][k] & DL[jpok][im1 - k]:
                    index_pair_1 = key >> shift
                    index_pair_2 = key & mask
                    parse_content = [index_pair_1, DPj[k][index_pair_1], index_pair_2, DP[jpok][im1 - k][index_pair_2]]
                    for rule_index in nt[key]:
                        DPji[rule_index] = parse_content
                        DL[j][i] |= lft[rule_index]
                        DR[j][i] |= rgt[rule_index]
                        if rule_index in substitution:
                            for s in substitution[rule_index]:
                                DPji[s] = parse_content
                                DL[j][i] |= lft[s]
                                DR[j][i] |= rgt[s]
                                    
            if DPji: Ksj.add(i)
            
            
    return DP