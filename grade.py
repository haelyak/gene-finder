import unittest
import find52 as hw

class SimpleTests(unittest.TestCase):
        
    def testB(self):
        # count
        shortTestSeq="ATGCCCCAGTGGTCAAGACGCCTCTCAACTGGGGAT"
        HWcodonCountD,HWtwoCodonCountD = hw.count([shortTestSeq])
        self.assertTrue(type(HWcodonCountD)==type({}) and type(HWtwoCodonCountD)==type({}),msg="count failed.")
        

    def testC(self):
        # count
        shortTestSeq="ATGCCCGAACCCCAA"
        HWcodonCountD,HWtwoCodonCountD  = hw.count([shortTestSeq])
        self.assertTrue(sum(HWtwoCodonCountD.values())==4,msg="count failed on "+ shortTestSeq)
        self.assertTrue(HWtwoCodonCountD["ATGCCC"]==1,msg="count failed on "+ shortTestSeq)
        self.assertTrue(HWtwoCodonCountD["CCCGAA"]==1,msg="count failed on "+ shortTestSeq)
        self.assertTrue(HWtwoCodonCountD["GAACCC"]==1,msg="count failed on "+ shortTestSeq)
        self.assertTrue(HWtwoCodonCountD["CCCCAA"]==1,msg="count failed on "+ shortTestSeq)

        self.assertTrue(sum(HWcodonCountD.values())==4,msg="count failed on "+ shortTestSeq)
        self.assertTrue(HWcodonCountD["ATG"]==1,msg="count failed on "+ shortTestSeq)
        self.assertTrue(HWcodonCountD["CCC"]==2,msg="count failed on "+ shortTestSeq)
        self.assertTrue(HWcodonCountD["GAA"]==1,msg="count failed on "+ shortTestSeq)

        
    def testD(self):
        # count
        HWcodonCountD,HWtwoCodonCountD = hw.count(hw.testSeqsL1)
        self.assertTrue(list(HWcodonCountD.values())==[61]*61,msg="count failed on testSeqsL1")
        self.assertTrue(list(HWtwoCodonCountD.values())==[1]*3721,msg="count failed on testSeqsL1")

        HWcodonCountD,HWtwoCodonCountD = hw.count(hw.testSeqsL2)
        self.assertTrue(sorted(HWcodonCountD.values())==54*[1]+5*[2]+[3,5],msg="count failed on testSeqsL2")
        self.assertTrue(sorted(HWtwoCodonCountD.values())==3649*[0]+72*[1],msg="count failed on testSeqsL2")

    def testE(self):
        # condProb
        HWprobD = hw.condProb(hw.testSeqsL1)
        self.assertTrue(type(HWprobD)==type({}))

        
    def testF(self):
        # condProb
        HWprobD = hw.condProb(hw.testSeqsL1)
        self.assertTrue([format(x,".3f") for x in HWprobD.values()]== ["0.016"]*3721,msg="condProb failed on testSeqsL1")  

    def testG(self):
        # condProb
        HWprobD=hw.condProb(hw.testSeqsL2)
        self.assertTrue([format(x,".2f") for x in sorted(HWprobD.values())]==["0.00"]*3649+["0.20"]*5+["0.33"]*3+["0.50"]*10+["1.00"]*54,msg="condProb failed on testSeqsL2")

    def testH(self):
        # logLikelihoodRatioSum
        codeProbD = hw.condProb(ecoliCodeOrfL)
        noncodeProbD = hw.condProb(ecoliNoncodeOrfL)
        llrD=hw.makeLogLikelihoodRatioD(codeProbD,noncodeProbD)
        self.assertTrue(type(hw.logLikelihoodRatioSum(hw.testOrf1,llrD))==type(1.1))

    def testI(self):
        # logLikelihoodRatioSum
        ecoliCodeProbD=hw.condProb(ecoliCodeOrfL)
        ecoliNoncodeProbD=hw.condProb(ecoliNoncodeOrfL)
        llrD=hw.makeLogLikelihoodRatioD(ecoliCodeProbD,ecoliNoncodeProbD)
        self.assertTrue(format(hw.logLikelihoodRatioSum(hw.testOrf1,llrD),".3f")==format(3.131517,".3f"),msg="logLikelihoodRatioSum failed on testOrf1")
        self.assertTrue(format(hw.logLikelihoodRatioSum(hw.testOrf2,llrD),".3f")==format(-22.340812,".3f"),msg="logLikelihoodRatioSum failed on testOrf2")

        
    def testJ(self):
        # predict
        codeProbD = hw.condProb(ecoliCodeOrfL)
        noncodeProbD = hw.condProb(ecoliNoncodeOrfL)
        llrD=hw.makeLogLikelihoodRatioD(codeProbD,noncodeProbD)
        predCodeL,predNonCodeL=hw.predict(vibAllOrfL[:100],llrD)
        self.assertTrue(len(predCodeL)==24 and len(predNonCodeL)==76,msg="predict failed on vibAllOrfL[:100]")


if __name__ == '__main__':
    ecoliCodeOrfL=[s.rstrip() for s in open("ecoliCodeTrainOrfs-mid.txt","r").readlines()]
    ecoliNoncodeOrfL=[s.rstrip() for s in open("ecoliNcTrainOrfs-mid.txt","r").readlines()]
    vibAllOrfL=[s.rstrip() for s in open("vibAllOrfs-mid.txt","r").readlines()]

    unittest.main()
