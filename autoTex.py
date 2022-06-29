from additiveCharacters import *
from dirichletCharacters import *
from dedekindSum import *
from SLTwoZ import *

def schreierGraphSLTwoZOverGammaZero(p, outFilePath, showST = [True, True]):
    ls = []
    ls.append(r'\begin{center}' + '\n')
    ls.append(r'\resizebox{\textwidth}{!}{'+ '\n')
    ls.append(r'\begin{tikzpicture}[node distance={20mm}, main/.style = {draw, circle}]'+ '\n')

    ls.append(r'\node[main] (1) {$1$};'+ '\n')
    ls.append(r'\node[main] (2) [right of=1] {$S$};'+ '\n')
    for i in range(1, p):
        ls.append(r'\node[main] (' + str(2 + i) + ') [right of=' + str (1 + i) + '] {$ST^{' + str(i) + '}$};'+ '\n')

    gaS = SLTwoZOverGammaZeroGroupAction(p, 'S')
    gaT = SLTwoZOverGammaZeroGroupAction(p, 'T')

    if showST[0]:
        visitedS = []
        ls.append(r'\draw [<->] (1) -- node[midway, above] {S} (2);'+ '\n')
        for i in range(1, p):
            if not (i in visitedS):
                ls.append(r'\draw [<->] (' + str(i + 2) + ') to [out=315, in=225, looseness=1]node[below]{S} (' + str(gaS[i] + 2) + ');'+ '\n')
                visitedS.append(i)
                visitedS.append(gaS[i])

    if showST[1]:
        visitedT = []
        ls.append(r'\draw [->] (1) edge[loop left]node{T} (1);'+ '\n')
        for i in range(0, p - 1):
            if not (i in visitedT):
                ls.append(r'\draw [->] (' + str(i + 2) + ') -- node[midway, above] {T} (' + str(gaT[i] + 2) + ');'+ '\n')
                visitedT.append(i)
        ls.append(r'\draw [<-] (2) to [out=300,in=240,looseness=1]node[below]{T} (' + str(p + 1) + ');'+ '\n')

    ls.append(r'\end{tikzpicture}'+ '\n')
    ls.append(r'}'+ '\n')
    ls.append(r'\end{center}'+ '\n')

    f = open(outFilePath, "w")
    f.writelines(ls)
    f.close()