from chimerax.core.commands import run

struc = {"HelixAlphaC":[728,729,730,731,732,733,734,735,736,737,738,739,740,741,742,743,744],
"HelixAlphaF":[868,867,868,869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885],
"RegulatorySpine":[753,742,832,811],
"CatalyticSpine":[702,819,820,821,774],
"ActivationLoop":[831,832,833,834,835,836,837,838,839,840,841,842,843,844,845,846,847,848,849,850,851,852],
"P-Loop":[695,696,697,698,699,700],
"Glu":[738],
"Lys":[721],
"DFG":[831],
"ASP":[813],
"ASP2":[870]
}

def intlisttostr(intlist):
    intstr=str(intlist)
    intstr=':'+intstr[1:len(intstr)-1]
    return intstr

def name_structures():
    for i in struc:
        run(session, 'name %s %s' %(i,intlisttostr(struc.get(i))))
def show_structures(structures):
    for i in struc:
        if i in structures:
            run(session, 'show %s' %(intlisttostr(struc.get(i))))
def style_structures(structures,style):
    for i in struc:
        if i in structures:
            run(session, 'style %s %s' %(intlisttostr(struc.get(i)),style))
def color_structures(structures,color,spec):
    for i in struc:
        if i in structures:
            run(session, 'color %s %s %s' %(intlisttostr(struc.get(i)),color,spec))
def transparency_structures(structures,percenttrans):
    for i in struc:
        if i in structures:
            run(session, 'transparency %s %s atoms' %(intlisttostr(struc.get(i)),percenttrans))
def label_structures(structures):
    for i in struc:
        if i in structures:
            run(session, 'label :%s text %s' %(struc.get(i)[int(len(struc.get(i))/2)],i))
def label_protein(session):
    run(session, 'color cyan ribbons')
    run(session, 'cartoon style modeHelix tube radius 1.5')
    name_structures()

    color_structures(['ActivationLoop'],'red','ribbon')
    color_structures(['P-Loop'],'pink','ribbons')

    show_structures(["Glu","Lys","DFG","ASP"])
    style_structures(["Glu","Lys","DFG","ASP"],"ball")
    color_structures(["Glu","Lys","DFG","ASP"],"grey",'atoms')
    color_structures(["Glu","Lys","DFG","ASP"],"byhetero",'atoms')

    show_structures(['RegulatorySpine','CatalyticSpine'])
    style_structures(['RegulatorySpine','CatalyticSpine'],"sphere")
    color_structures(['RegulatorySpine'],"purple",'atoms')
    color_structures(['CatalyticSpine'],"yellow",'atoms')
    transparency_structures(['RegulatorySpine','CatalyticSpine'],'50')

    label_structures(["Glu","Lys","DFG","ASP",'RegulatorySpine','CatalyticSpine','ActivationLoop','P-Loop','HelixAlphaC','HelixAlphaF'])
    run(session, 'view matrix camera 0.97361,0.17145,0.15062,40.948,0.1459,-0.97512,0.16687,44.571,0.17548,-0.14049,-0.97441,-163.76')


def register_command(session):
    from chimerax.core.commands import CmdDesc, register
    desc = CmdDesc(synopsis='Label EGFR protein')
    register('EGFRLabel', desc, label_protein, logger=session.logger)

register_command(session)
