import sympy as sym

# get the independent nonzero beta2 elements (sympy tensor chi2[xyz]/b2[abc]) from molecular symmetry
def nonzerob2(T,symmetry):
    tl = []
    if T == 'b2':
        ch = 'abc'
    else:
        ch = 'xyz'
    for i in ch:
        for j in ch:
            for k in ch:
                x = sym.symbols(T+i+j+k)
                tl.append(x)
    if symmetry == 'c1':
        tl = tl
    elif symmetry == 'c2':
        tl[0]=tl[4]=tl[8]=tl[6]=tl[2]=tl[14]=tl[16]=tl[10]=tl[12]=tl[18]=tl[22]=tl[26]=tl[24]=tl[20]=0
    elif symmetry == 'c3':
        tl[5]=tl[7]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[15]=tl[11]=tl[23]=tl[25]=tl[19]=tl[21]=0
    elif symmetry == 'd2':
        tl[0]=tl[4]=tl[8]=tl[6]=tl[2]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[14]=tl[16]=tl[10]=tl[12]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=0
    elif symmetry == 'c2v':
        tl[0]=tl[4]=tl[8]=tl[5]=tl[7]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[15]=tl[11]=tl[10]=tl[12]=tl[23]=tl[25]=tl[24]=tl[20]=tl[19]=tl[21]=0
    elif symmetry == 'c4':
        tl[0]=tl[4]=tl[8]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[10]=tl[12]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[14]=tl[2]
        tl[16]=tl[6]
        tl[15]=-tl[7]
        tl[11]=-tl[5]
        tl[22]=tl[18]
        tl[21]=-tl[19]
    elif symmetry == 's4':
        tl[0]=tl[4]=tl[8]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[10]=tl[12]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[14]=-tl[2]
        tl[16]=-tl[6]
        tl[15]=tl[7]
        tl[11]=tl[5]
        tl[22]=-tl[18]
        tl[21]=tl[19]
    elif symmetry == 'd4':
        tl[0]=tl[4]=tl[8]=tl[6]=tl[2]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[14]=tl[16]=tl[10]=tl[12]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[15]=-tl[7]
        tl[11]=-tl[5]
        tl[21]=-tl[19]
    elif symmetry == 'c4v':
        tl[0]=tl[4]=tl[8]=tl[5]=tl[7]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[15]=tl[11]=tl[10]=tl[12]=tl[23]=tl[25]=tl[24]=tl[20]=tl[19]=tl[21]=0
        tl[14]=tl[2]
        tl[16]=tl[6]
        tl[22]=tl[18]
    elif symmetry == 'd2d':
        tl[0]=tl[4]=tl[8]=tl[6]=tl[2]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[14]=tl[16]=tl[10]=tl[12]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[15]=tl[7]
        tl[11]=tl[5]
        tl[21]=tl[19]
    elif symmetry == 'c3':
        tl[8]=tl[17]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[4]=tl[10]=tl[12]=-tl[0]
        tl[1]=tl[3]=tl[9]=-tl[13]
        tl[11]=tl[5]
        tl[15]=-tl[7]
        tl[14]=tl[2]
        tl[16]=tl[6]
        tl[22]=tl[18]
        tl[21]=-tl[19]
    elif symmetry == 'd3':
        tl[8]=tl[6]=tl[2]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[14]=tl[16]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[4]=tl[10]=tl[12]=-tl[0]
        tl[15]=-tl[7]
        tl[11]=-tl[5]
        tl[21]=-tl[19]
    elif symmetry == 'c3v':
        tl[0]=tl[4]=tl[8]=tl[5]=tl[7]=tl[17]=tl[15]=tl[11]=tl[10]=tl[12]=tl[23]=tl[25]=tl[24]=tl[20]=tl[19]=tl[21]=0
        tl[1]=tl[3]=tl[9]=-tl[13]
        tl[14]=tl[2]
        tl[16]=tl[6]
        tl[22]=tl[18]
    elif symmetry == 'c6':
        tl[0]=tl[4]=tl[8]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[10]=tl[12]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[14]=tl[2]
        tl[16]=tl[6]
        tl[15]=-tl[7]
        tl[11]=-tl[5]
        tl[22]=tl[18]
        tl[21]=-tl[19]
    elif symmetry == 'c3h':
        tl[8]=tl[5]=tl[7]=tl[6]=tl[2]=tl[17]=tl[14]=tl[16]=tl[15]=tl[11]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=tl[19]=tl[21]=0
        tl[4]=tl[10]=tl[12]=-tl[0]
        tl[1]=tl[3]=tl[9]=-tl[13]
    elif symmetry == 'd6':
        tl[0]=tl[4]=tl[8]=tl[6]=tl[2]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[14]=tl[16]=tl[10]=tl[12]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[15]=-tl[7]
        tl[11]=-tl[5]
        tl[21]=-tl[19]
    elif symmetry == 'c6v':
        tl[0]=tl[4]=tl[8]=tl[5]=tl[7]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[15]=tl[11]=tl[10]=tl[12]=tl[23]=tl[25]=tl[24]=tl[20]=tl[19]=tl[21]=0
        tl[14]=tl[2]
        tl[16]=tl[6]
        tl[22]=tl[18]
    elif symmetry == 'd3h':
        tl[0]=tl[4]=tl[8]=tl[5]=tl[7]=tl[6]=tl[2]=tl[17]=tl[14]=tl[16]=tl[15]=tl[11]=tl[10]=tl[12]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=tl[19]=tl[21]=0
        tl[1]=tl[3]=tl[9]=-tl[13]
    elif symmetry == 'o':
        tl[0]=tl[4]=tl[8]=tl[6]=tl[2]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[14]=tl[16]=tl[10]=tl[12]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[15]=tl[19]=tl[5]
        tl[7]=tl[11]=tl[21]=-tl[5]
    elif symmetry == 'td':
        tl[0]=tl[4]=tl[8]=tl[6]=tl[2]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[14]=tl[16]=tl[10]=tl[12]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[15]=tl[19]=tl[7]=tl[11]=tl[21]=tl[5]
    elif symmetry == 't':
        tl[0]=tl[4]=tl[8]=tl[6]=tl[2]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[14]=tl[16]=tl[10]=tl[12]=tl[18]=tl[22]=tl[26]=tl[23]=tl[25]=tl[24]=tl[20]=0
        tl[15]=tl[19]=tl[5]
        tl[11]=tl[21]=tl[7]
    elif symmetry == 'civ':
        tl[0]=tl[4]=tl[8]=tl[5]=tl[7]=tl[1]=tl[3]=tl[9]=tl[13]=tl[17]=tl[15]=tl[11]=tl[10]=tl[12]=tl[23]=tl[25]=tl[24]=tl[20]=tl[19]=tl[21]=0
        tl[14]=tl[2]
        tl[16]=tl[6]
        tl[22]=tl[18]
    else:
        for i in range(len(tl)):
            tl[i]=0
    beta2 = sym.MutableDenseNDimArray(tl, (3,3,3))
    return beta2