
def parser_core(text_to_parse, dict nt, dict t):
    
    cdef:
        list dp_table, ld, Ks
        dict lld
        set s
        long rule_index, mask, new_key, old_key, i_shift, index_pair_1, index_pair_2, key
        int n = len(text_to_parse), i, j, k, jpok, im1
    dp_table = list()
    Ks = list()
    
    mask = (1 << 32) - 1
    for i in range(n):
        ld = []
        for j in range(n - i):
            lld = dict()
            ld.append(lld)
        dp_table.append(ld)
        
        
    for i in range(n):
        c = text_to_parse[i]
        if c not in t: return None
        
        for rule_index in t[c]:
            new_key = rule_index >> 32
            old_key = rule_index & mask
            dp_node = [c, None, old_key, None]
            dp_table[i][0][new_key] = dp_node
        s = set()
        s.add(0)
        Ks.append(s)
    
    
    for i in range (1, n):
        im1 = i - 1
        for j in range(n - i):
            
            D, jp1 = dp_table[j], j + 1
            Di = D[i]
            Ksj = Ks[j]
            
            for k in Ksj:
                jpok = jp1 + k
                if im1 - k in Ks[jpok]:                
                    D1, D2 = D[k], dp_table[jpok][im1 - k]
                
                    for index_pair_1 in D1:
                        
                        i_shift = index_pair_1 * 4294967296
                        for index_pair_2 in D2:
                            key = i_shift + index_pair_2
                            
                            if key in nt:
                                parse_content = [index_pair_1, D1[index_pair_1], index_pair_2, D2[index_pair_2]]
                                for rule_index in nt[key]:
                                    Di[rule_index] = parse_content
            
            if len(D[i]) > 0: Ks[j].add(i)
            
            
    return dp_table