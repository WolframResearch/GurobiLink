
Needs["GurobiLink`"]

VerificationTest[
	GurobiTestLicense[]
	,
	True
	,
	TestID->"GurobiLink_Basic-20210519-C6Z5H7"
]

VerificationTest[
	(data = GurobiDataCreate[]) // Head
	,
	GurobiData
	,
	TestID->"GurobiLink_Basic-20210519-C3D2B1"
]

VerificationTest[
	GurobiDataQ[data]
	,
	True
	,
	TestID->"GurobiLink_Basic-20210519-I3N7Q8"
]

VerificationTest[
	GurobiSetVariableTypesAndObjectiveVector[data, {}, {1.}]
	,
	0
	,
	TestID->"GurobiLink_Basic-20210519-R4F6X5"
]

VerificationTest[
	GurobiAddLinearConstraints[data, SparseArray[{{1}}], ">", {0.}]
	,
	0
	,
	TestID->"GurobiLink_Basic-20210519-F6I2F8"
]

VerificationTest[
	GurobiOptimize[data]
	,
	"Solved"
	,
	TestID->"GurobiLink_Basic-20210519-W0L9Y9"
]

VerificationTest[
	GurobiObjectiveValue[data]
	,
	0.
	,
	TestID->"GurobiLink_Basic-20210519-D6R9J9"
]

VerificationTest[
	Gurobix[data]
	,
	{0.}
	,
	TestID->"GurobiLink_Basic-20210519-U2T2Z1"
]
