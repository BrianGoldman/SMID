#!/usr/bin/env python
from openeye.oechem import *
from openeye.oeomega import OEOmega, OEFlipper
import itertools
import math
import logging
import sys
from docopt import docopt

""" Set up the logging """
logger = logging.getLogger()
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s: %(levelname)s : %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


class SmidCalculator:
    def __init__(self, maxConfs=200, energyWindow=5., strictStereo=False):
        self.omega = self._getOmega(maxConfs, energyWindow, strictStereo)
        logger.info("Createing confs with maxConfs = %d eWindow = %f ", maxConfs, energyWindow)

    def calcSmid(self, inpMol, log=True):
        if not self.omega(inpMol):
            logger.info("Omega Failed on molecule %s", mol.GetTitle())
            raise Exception("Omega Fail")
        if log:
            logger.info("Created %d confs and now running SMID Calc ", inpMol.NumConfs())
        smidForMol = min(map(self._calcSmidForConf, inpMol.GetConfs()))
        return math.sqrt(smidForMol)

    def smidHScan(self, inpMol, verbose=True):
        mol = OEGraphMol(inpMol)
        seenSet = set()
        res = []

        def getNbr(a):
            for nay in a.GetAtoms():
                return nay

        hAtms = [_ for _ in mol.GetAtoms() if _.GetAtomicNum() == 1]
        if verbose:
            logger.info("Staring hScan on %d positions ", len(hAtms))

        for a in hAtms:
            # a is the atom on the original molecule 
            # find its neighbor and store the idx
            attchIdx = 0
            for attchAtm in a.GetAtoms():
                attchIdx = attchAtm.GetIdx()

            newMol = OEGraphMol(mol)
            tmpMol = OEGraphMol()
            OESmilesToMol(tmpMol, "c1c[c:1]ccc1")
            OEAddMols(newMol, tmpMol)
            atms = [_ for _ in newMol.GetAtoms() if _.GetMapIdx() == 1 or _.GetIdx() == a.GetIdx()]
            newMol.NewBond(getNbr(atms[0]), atms[1], 1)
            hIdx = atms[0].GetIdx()
            newMol.DeleteAtom(atms[0])
            atms[1].SetMapIdx(0)

            smi = OEMolToSmiles(newMol)
            if smi not in seenSet:
                seenSet.add(smi)
                n = len(list(OEFlipper(newMol, 2, False, False)))
                for enat in OEFlipper(newMol, 2, False, False):
                    if verbose and n > 1:
                        logger.info("Doing Enatomeric Molecule %s ", OEMolToSmiles(enat))
                    hSmid = self.calcSmid(OEMol(enat), log=False)
                    res.append((hIdx, attchIdx, hSmid))
        return res

    @staticmethod
    def _atmDist(atmPair):
        (c1, c2) = atmPair
        return math.pow(c1[0] - c2[0], 2) + math.pow(c1[1] - c2[1], 2) + math.pow(c1[2] - c2[2], 2)

    @staticmethod
    def _calcSmidForConf(conf):
        return max(map(SmidCalculator._atmDist, itertools.combinations(conf.GetCoords().values(), 2)))

    @staticmethod
    def _getOmega(maxConfs, energyWindow, strictStereo):
        omega = OEOmega()
        omega.SetMaxConfs(maxConfs)
        omega.SetFromCT(True)
        omega.SetStrictStereo(strictStereo)
        omega.SetEnergyWindow(energyWindow)
        return omega


if __name__ == "__main__":
    __doc__ = """smid.py
       Usage:
        smid.py -i=<file>  [ -o=<file> -n=<int> -e=<float> -s ]

        -i --input=<file>       Input file
        -o --output=<file>      Output file
        -n --numConfs=<int>     Number of confs [default: 200]
        -e --eWindow=<float>    Energy Window for omega [default: 5.0]
        -s, --scan              Do the HScan with attached benzene ring
    """

    arguments = docopt(__doc__)
    inp = oemolistream(arguments["--input"])
    fOut = sys.stdout if arguments["--output"] is None else open(arguments["--output"], "w")
    if arguments["--scan"]:
        fOut.write("Name,smid,h_idx,attach_idx,hsmid\n")
    else:
        fOut.write("Name,smid\n")

    smidCalculator = SmidCalculator(maxConfs=int(arguments["--numConfs"]), energyWindow=float(arguments["--eWindow"]))
    for mol in inp.GetOEMols():
        OEAddExplicitHydrogens(mol)
        # noinspection PyBroadException
        try:
            smid = smidCalculator.calcSmid(OEMol(mol))
            if arguments["--scan"]:
                hRes = smidCalculator.smidHScan(mol)
                fOut.write(f"{mol.GetTitle()},{smid:5.2f},,,\n")
                for h_idx, attach_idx, h_smid in hRes:
                    fOut.write(f"{mol.GetTitle()},,{h_idx},{attach_idx},{h_smid}\n")
            else:
                fOut.write(f"{mol.GetTitle()},{smid:5.2f}\n")
        except Exception as e:
            logger.info(f"Skipping molecule {mol.GetTitle()}")
            logger.info(f"{e}")
