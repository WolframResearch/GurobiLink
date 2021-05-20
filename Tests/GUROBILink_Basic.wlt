
Needs["GUROBILink`"]

VerificationTest[
	GUROBITestLicense[]
	,
	True
	,
	TestID->"GUROBILink_Basic-20210519-C6Z5H7"
]

VerificationTest[
	(data = GUROBIDataCreate[]) // Head
	,
	GUROBIData
	,
	TestID->"GUROBILink_Basic-20210519-C3D2B1"
]

VerificationTest[
	GUROBIDataQ[data]
	,
	True
	,
	TestID->"GUROBILink_Basic-20210519-I3N7Q8"
]

VerificationTest[
	GUROBISetVariableTypesAndObjectiveVector[data, {}, {1.}]
	,
	0
	,
	TestID->"GUROBILink_Basic-20210519-R4F6X5"
]

VerificationTest[
	GUROBIAddLinearConstraints[data, SparseArray[{{1}}], ">", {0.}]
	,
	0
	,
	TestID->"GUROBILink_Basic-20210519-F6I2F8"
]

VerificationTest[
	GUROBIOptimize[data]
	,
	"Solved"
	,
	TestID->"GUROBILink_Basic-20210519-W0L9Y9"
]

VerificationTest[
	GUROBIObjectiveValue[data]
	,
	0.
	,
	TestID->"GUROBILink_Basic-20210519-D6R9J9"
]

VerificationTest[
	GUROBIx[data]
	,
	{0.}
	,
	TestID->"GUROBILink_Basic-20210519-U2T2Z1"
]
