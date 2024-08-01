


'''

用于生成v0,v1和nonterminal_nodes的新图
mapping:新点的映射


输入:
    *k: 最短路径集合的节点集
    *nonpathnodes:  不在k集合中的节点

* mapping:  m0,m1,m2,v3,v7,v9,v15
            0,1,2,3,4,5,6

'''

def st_3(k,nonpathnodes):

    nodes_set = k + nonpathnodes
    l = len(k) + len(nonpathnodes)

    new_nodes = [i for i in range(l)]
    mapping_nodes =  zip([i for i in range(l)],nodes_set)
    mapping = [list(x) for x in mapping_nodes]

    mapping = dict(mapping)

    #print('st_3:将P_diff的节点(m0,m1..)和剩余节点映射到新图上  mapping:',mapping)

    return mapping