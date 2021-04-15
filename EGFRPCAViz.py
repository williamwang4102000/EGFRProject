def PCA_Dist12(session):
    from chimerax.core.commands import all_objects
    m=session.models[0]
    clust1=[[705,778],[767,776],[745,855]]
    clust2=[[791,852],[764,768],[766,790],[759,788],[762,788],[762,797],[745,759],[745,762]]
    for i in clust1:
        residue1=m.find_residue('A',i[0])
        residue2=m.find_residue('A',i[1])
        get_shortest_distance(session, residue1, residue2, 'green')
    for j in clust2:
        residue1=m.find_residue('A',j[0])
        residue2=m.find_residue('A',j[1])
        get_shortest_distance(session, residue1, residue2, 'red')
def PCA_Dist23(session):
    from chimerax.core.commands import all_objects
    m=session.models[0]
    clust2=[[768,831],[758,857],[765,856],[853,856],[774,856],[837,855],[842,855]]
    clust3=[[775,855],[792,844],[790,855],[743,856],[793,794]]
    for i in clust2:
        residue1=m.find_residue('A',i[0])
        residue2=m.find_residue('A',i[1])
        get_shortest_distance(session, residue1, residue2, 'red')
    for j in clust3:
        residue1=m.find_residue('A',j[0])
        residue2=m.find_residue('A',j[1])
        get_shortest_distance(session, residue1, residue2, 'cyan')
def PCA_Dist24(session):
    from chimerax.core.commands import all_objects
    m=session.models[0]
    clust2=[[774,856],[759,788],[764,768],[828,833],[765,856],[768,831],[791,852]]
    clust4=[[771,773],[740,741],[745,855]]
    for i in clust2:
        residue1=m.find_residue('A',i[0])
        residue2=m.find_residue('A',i[1])
        get_shortest_distance(session, residue1, residue2, 'red')
    for j in clust4:
        residue1=m.find_residue('A',j[0])
        residue2=m.find_residue('A',j[1])
        get_shortest_distance(session, residue1, residue2, 'purple')
def get_shortest_distance(session, residue1, residue2, color):
    import math
    atoms1=residue1.atoms
    atoms2=residue2.atoms
    min_dis=1000000000000.0
    min_atoms_1 = atoms1[0]
    min_atoms_2= atoms2[0]
    for a in atoms1:
        for b in atoms2:
            dis =math.sqrt((a.coord[0]-b.coord[0])**2+(a.coord[1]-b.coord[1])**2+(a.coord[2]-b.coord[2])**2)
            if(dis<min_dis):
                min_dis=dis
                min_atoms_1=a
                min_atoms_2=b
    session.logger.info(str(min_atoms_1))
    session.logger.info(str(min_atoms_2))
    session.logger.info(str(min_dis))
    sma1=str(min_atoms_1).split()
    sma2=str(min_atoms_2).split()
    session.logger.info(str(sma1))
    from chimerax.core.commands import run
    run(session, 'show :' + sma1[2]+' atoms')
    run(session, 'show :' + sma2[2]+' atoms')
    run(session, 'distance '+sma1[0]+':'+sma1[2]+'@'+sma1[3]+' '+sma2[0]+':'+sma2[2]+'@'+sma2[3]+' color '+color)

def register_command(session):
    from chimerax.core.commands import CmdDesc, register, StringArg
    desc12 = CmdDesc(synopsis='vizualize egfr12 distances')
    desc23 = CmdDesc(synopsis='vizualize egfr23 distances')
    desc24 = CmdDesc(synopsis='vizualize egfr24 distances')
    register('PCADist12', desc12, PCA_Dist12, logger=session.logger)
    register('PCADist23', desc23, PCA_Dist23, logger=session.logger)
    register('PCADist24', desc24, PCA_Dist24, logger=session.logger)
register_command(session)
