
Needs["GUROBILink`"]

Test[
	GUROBITestLicense[]
	,
	True
	,
	TestID->"GUROBILink_Basic-20210519-C6Z5H7"
]

Test[
	(data = GUROBIDataCreate[]) // Head
	,
	GUROBIData
	,
	TestID->"GUROBILink_Basic-20210519-C3D2B1"
]

Test[
	GUROBIDataQ[data]
	,
	True
	,
	TestID->"GUROBILink_Basic-20210519-I3N7Q8"
]

Test[
	GUROBISetVariableTypesAndObjectiveVector[data, {}, {1.}]
	,
	0
	,
	TestID->"GUROBILink_Basic-20210519-R4F6X5"
]

Test[
	GUROBIAddLinearConstraints[data, SparseArray[{{1}}], ">", {0.}]
	,
	0
	,
	TestID->"GUROBILink_Basic-20210519-F6I2F8"
]

Test[
	GUROBIOptimize[data]
	,
	"Solved"
	,
	TestID->"GUROBILink_Basic-20210519-W0L9Y9"
]

Test[
	GUROBIObjectiveValue[data]
	,
	0.
	,
	TestID->"GUROBILink_Basic-20210519-D6R9J9"
]

Test[
	GUROBIx[data]
	,
	{0.}
	,
	TestID->"GUROBILink_Basic-20210519-U2T2Z1"
]
