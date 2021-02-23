# Create command to save an image of each chain of an atomic structure.
#
#  chainimages #1
#
def label_protein(session, protein):
    from chimerax.core.commands import run
    if(protein.upper()=='1ATP'):
        #opens protein model
        run(session, 'open '+protein)

        #Names key structures
        run(session, 'name HelixAlphac :85-97')#
        run(session, 'name Glu91 :91')#
        run(session, 'name ActivationLoop	:184-201')#
        run(session, 'name D166 :166')#
        run(session, 'name P-Loop :51-55')#
        run(session, 'name Lys72 :72')
        run(session, 'name RegulatorySpine	:106,95,185,164')#
        run(session, 'name CatalyticSpine :231, 227, 128, 172, 174, 173, 57, 70')#
        run(session, 'name p-T197 :197')#
        run(session, 'name DFG184 :184')
        run(session, 'name D166 :166')
        run(session, 'name HelixAlphaf :218-233')
        run(session, 'name ASP220 :220')

    #show or hide sidechains/atoms of structures
        run(session, 'show Glu91')
        run(session, 'show p-T197')
        run(session, 'show D166')
        run(session, 'hide P-Loop')
        run(session, 'show RegulatorySpine')
        run(session, 'show CatalyticSpine')
        run(session, 'show ASP220')
    #style structures
        run(session, 'cartoon style modeHelix tube sides 10')
        run(session, 'style ligand ball')
        run(session, 'style Glu91 ball')
        run(session, 'style p-T197 ball')
        run(session, 'style D166 ball')
        run(session, 'style RegulatorySpine sphere')
        run(session, 'style CatalyticSpine sphere')
        run(session, 'style Lys72 ball')
        run(session, 'style DFG184 ball')
        run(session, 'style D166 ball')
        run(session, 'style ASP220')
    #color structures
        run(session, 'color cyan')
        run(session, 'color HelixAlphac blue cartoons')
        run(session, 'color HelixAlphaf lightblue cartoons')
        run(session, 'color ::name="MN" lime')
        run(session, 'color ligand orange')
        run(session, 'color ligand byhetero')
        run(session, 'color ActivationLoop red')
        run(session, 'color P-Loop pink')

        run(session, 'color p-T197 grey atoms')
        run(session, 'color p-T197 byhetero atoms')
        run(session, 'color Glu91 grey atoms')
        run(session, 'color Glu91 byhetero atoms')
        run(session, 'color D166 grey atoms')
        run(session, 'color D166 byhetero atoms')
        run(session, 'color RegulatorySpine purple atoms')
        run(session, 'transparency RegulatorySpine 50 atoms')
        run(session, 'color CatalyticSpine yellow atoms')
        run(session, 'transparency CatalyticSpine 50 atoms')
        run(session, 'color Lys72 grey atoms')
        run(session, 'color Lys72 byhetero atoms')
        run(session, 'color DFG184 grey atoms')
        run(session, 'color DFG184 byhetero atoms')
        run(session, 'color D166 grey atoms')
        run(session, 'color D166 byhetero atoms')
        run(session, 'color ASP220 grey atoms')
        run(session, 'color ASP220 byhetero atoms')
    elif(protein.upper()=='1M14'):
        run(session, 'open 1M14')
        #Names key structures
        run(session, 'name HelixAlphac :728-744')#
        run(session, 'name Glu91 :738')#
        run(session, 'name ActivationLoop	:831-852')#
        run(session, 'name D166 :813')#
        run(session, 'name P-Loop :695-700')#
        run(session, 'name Lys72 :721')
        run(session, 'name RegulatorySpine	:753,742,832,811')#
        run(session, 'name CatalyticSpine :726, 798, 843, 844')#
        run(session, 'name DFG184 :831')
        run(session, 'name D166 :813')
        run(session, 'name HelixAlphaf :868-885')
        run(session, 'name ASP220 :870')

    #show or hide sidechains/atoms of structures
        run(session, 'show Glu91')
        run(session, 'show D166')
        run(session, 'hide P-Loop')
        run(session, 'show RegulatorySpine')
        run(session, 'show CatalyticSpine')
        run(session, 'show ASP220')
    #style structures
        run(session, 'cartoon style modeHelix tube sides 10')
        run(session, 'style Glu91 ball')
        run(session, 'style D166 ball')
        run(session, 'style RegulatorySpine sphere')
        run(session, 'style CatalyticSpine sphere')
        run(session, 'style Lys72 ball')
        run(session, 'style DFG184 ball')
        run(session, 'style D166 ball')
        run(session, 'style ASP220')
    #color structures
        run(session, 'color cyan')
        run(session, 'color HelixAlphac blue cartoons')
        run(session, 'color HelixAlphaf lightblue cartoons')
        run(session, 'color ActivationLoop red')
        run(session, 'color P-Loop pink')

        run(session, 'color Glu91 grey atoms')
        run(session, 'color Glu91 byhetero atoms')
        run(session, 'color D166 grey atoms')
        run(session, 'color D166 byhetero atoms')
        run(session, 'color RegulatorySpine purple atoms')
        run(session, 'transparency RegulatorySpine 50 atoms')
        run(session, 'color CatalyticSpine yellow atoms')
        run(session, 'transparency CatalyticSpine 50 atoms')
        run(session, 'color Lys72 grey atoms')
        run(session, 'color Lys72 byhetero atoms')
        run(session, 'color DFG184 grey atoms')
        run(session, 'color DFG184 byhetero atoms')
        run(session, 'color D166 grey atoms')
        run(session, 'color D166 byhetero atoms')
        run(session, 'color ASP220 grey atoms')
        run(session, 'color ASP220 byhetero atoms')

def register_command(session):
    from chimerax.core.commands import CmdDesc, register
    from chimerax.core.commands import StringArg
    desc = CmdDesc(required = [('protein', StringArg)], synopsis='Label EGFR protein')
    register('KinaseLabel', desc, label_protein, logger=session.logger)

register_command(session)
