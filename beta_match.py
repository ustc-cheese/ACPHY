from beamline import parseLattice, ParseParams, BeamMatch, FELSimulator

def test():
    testParse = ParseParams('rad.in', 'rad.lat')
    aw0 = testParse.getUndulatorParameter()
    xlamd = testParse.getUndulatorPeriod()
    unitlength = testParse.getUndulatorUnitlength()
    xlamds = testParse.getFELwavelength()
    gamma0 = testParse.getElectronGamma()
    emitx = testParse.getElectronEmitx()
    # imagl = testParse.getChicaneMagnetLength()
    # idril = testParse.getChicaneDriftLength()
    # ibfield = testParse.getChicaneMagnetField()

    print("aw0    = %.3f" % aw0   )
    print("xlamd  = %.3f" % xlamd )
    print("xlamds = %.3e" % xlamds)
    print("gamma0 = %.3f" % gamma0)
    print("emitx  = %.3e" % emitx )
    # print("imagl  = %.3f" % imagl )
    # print("idril  = %.3f" % idril )
    # print("ibfield= %.3f" % ibfield)
    print("unit   = %.3f" % unitlength)
    # print parseLattice('fullat.hghg')

    qf, qd = -1, 2

    testMatch = BeamMatch('mod.in', 'rad.in', 'mod.lat', 'rad.lat',\
                          'newmod.in', 'newrad.in', 'newrad.lat', qf, qd)
    if testMatch.matchCalculate():
        testMatch.matchPerform()
        testMatch.matchPrintout()
        fel = FELSimulator()
        fel.run()
        fel.postProcess()
        print(fel.getMaxPower())

if __name__ == '__main__':
        test()
