PROGRAM index
IMPLICIT NONE
CHARACTER(16) :: Molecule, Residue, ResidueType, AtomType
CHARACTER(128) :: line
INTEGER(kind=16) :: limit
INTEGER :: i, j, nLines, nHeader, ios
INTEGER :: AtomIndex, ResidueIndex, nResidues
INTEGER, ALLOCATABLE :: IndexStore(:,:)

! Read input values
	OPEN(0,file='./system.inp')
	READ(0,*)Molecule
	READ(0,*)nResidues
	CLOSE(0)

! Allocate IndexStore
	ALLOCATE(IndexStore(nResidues,8))

! Determine number of lines in topology file
	limit=1e10
	nLines=0
	ios = 0
	OPEN(1,file=TRIM(ADJUSTL(Molecule))//".gro")
	DO i=1,limit
		READ(1,*,IOSTAT=ios)line
		IF(ios /=0) EXIT
		IF (J == limit) THEN
			WRITE(*,*)"MAXLINES EXCEEDED"
			STOP
		END IF
		nLines=nLines+1
	END DO
	CLOSE(1)

	WRITE(*,*)
	WRITE(*,'(2A, I10, A)')TRIM(Molecule), ".gro has", nLines, " lines."
	WRITE(*,'(2A, I5, A)')TRIM(Molecule), " has ", nResidues, " residues."
	WRITE(*,*)

! Read atom index from topology file
	OPEN(10,file=TRIM(Molecule)//".gro")

! Skip topology header
! NOTE: This is designed for *.gro format which have 2 lines. It may be
! necessary to change for varying formats or from different format
! implementations.
	nHeader = 2
	DO i=1,nHeader
		READ(10,*)
	END DO

! Read topology information
	DO i=1,nLines-(nHeader+1)
		READ(10,'(A)')line
		READ(line(6:8),'(A)')ResidueType
		READ(line(12:16),'(A)')AtomType
		READ(line(16:20),'(I5)')AtomIndex

! Adenosine Residues
		IF(ResidueType=="DA ".OR.ResidueType=="DA5".OR.ResidueType=="DA3")THEN
			IF(i == 1)THEN
				ResidueIndex = 1
				Residue = ResidueType
				WRITE(*,'(A)')"Adenosine"
			ELSEIF(Residue/=ResidueType)THEN
				ResidueIndex = ResidueIndex + 1
				Residue = ResidueType
				WRITE(*,'(A)')"Adenosine"
			END IF

			IF(AtomType=="  C5 ")THEN
				IndexStore(ResidueIndex, 1) = AtomIndex
			ELSEIF(AtomType=="  C4 ")THEN
				IndexStore(ResidueIndex, 2) = AtomIndex
				IndexStore(ResidueIndex, 8) = AtomIndex
			ELSEIF(AtomType=="  N3 ")THEN
				IndexStore(ResidueIndex, 3) = AtomIndex
			ELSEIF(AtomType=="  N9 ")THEN
				IndexStore(ResidueIndex, 4) = AtomIndex
				IndexStore(ResidueIndex, 7) = AtomIndex
			ELSEIF(AtomType==" O4' ")THEN
				IndexStore(ResidueIndex, 5) = AtomIndex
			ELSEIF(AtomType==" C1' ")THEN
				IndexStore(ResidueIndex, 6) = AtomIndex
			END IF

! Cytidine Residues
		ELSEIF(ResidueType=="DC ".OR.ResidueType=="DC5".OR.ResidueType=="DC3")THEN
			IF(i == 1)THEN
				ResidueIndex = 1
				Residue = ResidueType
				WRITE(*,'(A)')"Cytosine"
			ELSEIF(Residue/=ResidueType)THEN
				ResidueIndex = ResidueIndex + 1
				Residue = ResidueType
				WRITE(*,'(A)')"Cytosine"
			END IF

			IF(AtomType=="  C4 ")THEN
				IndexStore(ResidueIndex, 1) = AtomIndex
			ELSEIF(AtomType=="  N1 ")THEN
				IndexStore(ResidueIndex, 2) = AtomIndex
				IndexStore(ResidueIndex, 7) = AtomIndex
			ELSEIF(AtomType=="  C2 ")THEN
				IndexStore(ResidueIndex, 3) = AtomIndex
				IndexStore(ResidueIndex, 8) = AtomIndex
			ELSEIF(AtomType=="  C6 ")THEN
				IndexStore(ResidueIndex, 4) = AtomIndex
			ELSEIF(AtomType==" C1' ")THEN
				IndexStore(ResidueIndex, 5) = AtomIndex
			ELSEIF(AtomType==" O4' ")THEN
				IndexStore(ResidueIndex, 6) = AtomIndex
			END IF

! Guanosine Residues
		ELSEIF(ResidueType=="DG ".OR.ResidueType=="DG5".OR.ResidueType=="DG3")THEN
			IF(i == 1)THEN
				ResidueIndex = 1
				Residue = ResidueType
				WRITE(*,'(A)')"Guanosine"
			ELSEIF(Residue/=ResidueType)THEN
				ResidueIndex = ResidueIndex + 1
				Residue = ResidueType
				WRITE(*,'(A)')"Guanosine"
			END IF

			IF(AtomType=="  C5 ")THEN
				IndexStore(ResidueIndex, 1) = AtomIndex
			ELSEIF(AtomType=="  C4 ")THEN
				IndexStore(ResidueIndex, 2) = AtomIndex
				IndexStore(ResidueIndex, 8) = AtomIndex
			ELSEIF(AtomType=="  N3 ")THEN
				IndexStore(ResidueIndex, 3) = AtomIndex
			ELSEIF(AtomType=="  N9 ")THEN
				IndexStore(ResidueIndex, 4) = AtomIndex
				IndexStore(ResidueIndex, 7) = AtomIndex
			ELSEIF(AtomType==" O4' ")THEN
				IndexStore(ResidueIndex, 5) = AtomIndex
			ELSEIF(AtomType==" C1' ")THEN
				IndexStore(ResidueIndex, 6) = AtomIndex
			END IF

! Thymidine Residues
		ELSE IF(ResidueType=="DT ".OR.ResidueType=="DT5".OR.ResidueType=="DT3")THEN
			IF(i == 1)THEN
				ResidueIndex = 1
				Residue = ResidueType
				WRITE(*,'(A)')"Thymidine"
			ELSEIF(Residue/=ResidueType)THEN
				ResidueIndex = ResidueIndex + 1
				Residue = ResidueType
				WRITE(*,'(A)')"Thymidine"
			END IF

			IF(AtomType=="  C4 ")THEN
				IndexStore(ResidueIndex, 1) = AtomIndex
			ELSEIF(AtomType=="  N1 ")THEN
				IndexStore(ResidueIndex, 2) = AtomIndex
				IndexStore(ResidueIndex, 7) = AtomIndex
			ELSEIF(AtomType=="  C2 ")THEN
				IndexStore(ResidueIndex, 3) = AtomIndex
				IndexStore(ResidueIndex, 8) = AtomIndex
			ELSEIF(AtomType=="  C6 ")THEN
				IndexStore(ResidueIndex, 4) = AtomIndex
			ELSEIF(AtomType==" C1' ")THEN
				IndexStore(ResidueIndex, 5) = AtomIndex
			ELSEIF(AtomType==" O4' ")THEN
				IndexStore(ResidueIndex, 6) = AtomIndex
			END IF
		END IF
	END DO

! Create atomic index files
	OPEN(100,file="./Atoms.ndx")
	WRITE(100,'(A)')"[ Atomic Sites ]"

! Write atomic index to index file
	DO j=1,nResidues
		DO i=1,8
			WRITE(100,*)IndexStore(j,i)
		END DO
	END DO
	CLOSE(100)

! Deallocate IndexStore
	DEALLOCATE(IndexStore)

END PROGRAM index
