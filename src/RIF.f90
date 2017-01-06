!###############################################################################
!                                                                              #
!         CCCCCC  SSSSSS  IIIII  RRRRRR   OOOO          L      IIIII           #
!         C       S         I    R    R  O    O         L        I             #
!         C       SSSSSS    I    RRRRRR  O    O  =====  L        I             #
!         C            S    I    R RR    O    O         L        I             #
!         CCCCCC  SSSSSS  IIIII  R   RR   OOOO          LLLLL  IIIII           #
!                                                                              #
!         Queensland Bioscience Precint                                        #
!         306 Carmody Road                                                     #
!         St Lucia, QLD 4067                                                   #
!         AUSTRALIA                                                            #
!         Ph. (07) 3214 2392                                                   #
!         Fx. (07) 3214 2900                                                   #
!                                                                              #
!         Objective:      Perform RIF algorithm                                #
!         Authors:        Antonio Reverter (Tony.Reverter-Gomez@csiro.au)      #
!                         Andrea Rau (andrea.rau@inra.fr)                      #
!						  Florence Jaffrezic (florence.jaffrezic@inra.fr)      #
!         Last Revision:  January 2016                                         #
!                                                                              #
!###############################################################################

SUBROUTINE rif(ntf, nta, nconditions1, nconditions2, &
	TFdata1, targetdata1, TFdata2, targetdata2, rif1, rif2, &
	cond1Xntf, cond2Xntf, cond1Xnta, cond2Xnta)

	IMPLICIT NONE
	! Type declarations
	INTEGER :: nconditions1, nconditions2, ntf, nta
	INTEGER :: cond1Xntf
	INTEGER :: cond2Xntf
	INTEGER :: cond1Xnta
	INTEGER :: cond2Xnta
	
	DOUBLE PRECISION :: TFdata1(cond1Xntf), TFdata2(cond2Xntf), targetdata1(cond1Xnta), targetdata2(cond2Xnta)
	DOUBLE PRECISION :: rif1(ntf), rif2(ntf)
	DOUBLE PRECISION :: tfcexpr(ntf, nconditions1), tfnexpr(ntf, nconditions2), tacexpr(nta, nconditions1), tanexpr(nta, nconditions2)
	DOUBLE PRECISION :: gene_ccorr(ntf, nta), gene_ncorr(ntf, nta)
	DOUBLE PRECISION :: tf_tmp1(nconditions1), ta_tmp1(nconditions1), tf_tmp2(nconditions2), ta_tmp2(nconditions2)

	DOUBLE PRECISION :: corr, ave, de, dw, er1, er2
	DOUBLE PRECISION :: tmp1, tmp2, tmp3, tmp4
	INTEGER :: i,j,k,ind,low1,high1,low2,high2
	
	!#################################################################
	

	DO i = 1,ntf
		low1=(i-1)*nconditions1 + 1
		high1=(i*nconditions1)
		low2=(i-1)*nconditions2 + 1
		high2=(i*nconditions2)
		ind = 1
		DO j = low1,high1
			tfcexpr(i,ind) = TFdata1(j)
			ind = ind + 1
		ENDDO
		ind = 1
		DO j = low2,high2
			tfnexpr(i,ind) = TFdata2(j)
			ind = ind + 1
		ENDDO
	ENDDO
	
	DO i = 1,nta
		low1=(i-1)*nconditions1 + 1
		high1=(i*nconditions1)
		low2=(i-1)*nconditions2 + 1
		high2=(i*nconditions2)
		ind = 1
		DO j = low1, high1
			tacexpr(i,ind) = targetdata1(j)
			ind = ind + 1
		ENDDO
		ind = 1
		DO j = low2, high2
			tanexpr(i,ind) = targetdata2(j)
			ind = ind + 1
		ENDDO
	ENDDO	
	
	!###################################################
	! Double loop and printing output
	!###################################################
		
	DO i = 1, ntf
		rif1(i) = 0.
		rif2(i) = 0.
		DO j = 1, nta
			tf_tmp1(:) = tfcexpr(i,:)
			ta_tmp1(:) = tacexpr(j,:)
			gene_ccorr(i,j) = corr(tf_tmp1, ta_tmp1, nconditions1)
			tf_tmp2(:) = tfnexpr(i,:)
			ta_tmp2(:) = tanexpr(j,:)
			gene_ncorr(i,j) = corr(tf_tmp2, ta_tmp2, nconditions2)
			ave = (SUM(ta_tmp1)/nconditions1 + SUM(ta_tmp2)/nconditions2)/2.
			de = SUM(ta_tmp1)/nconditions1 - SUM(ta_tmp2)/nconditions2
			dw = gene_ccorr(i,j) - gene_ncorr(i,j)
			rif1(i) = rif1(i) + ave * de * (dw**2)
			
			er1= SUM(ta_tmp1)/nconditions1 * gene_ccorr(i,j)
			er2= SUM(ta_tmp2)/nconditions2 * gene_ncorr(i,j)
			rif2(i) = rif2(i) + er1**2 - er2**2 
		ENDDO
		rif1(i) = rif1(i) / nta
		rif2(i) = rif2(i) / nta
	ENDDO
	
END SUBROUTINE 


!#########################################################################
!#########################################################################
FUNCTION corr(x,y,n) RESULT(f_val)
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: n
	DOUBLE PRECISION, INTENT(IN) :: x(n),y(n)
	DOUBLE PRECISION :: f_val
	DOUBLE PRECISION :: s1, ss1, s2, ss2, ss12, var1, var2, num, den
	INTEGER :: i,j

	s1=0.0; ss1=0.0; s2=0.0; ss2=0.0; ss12=0.0

	DO i = 1, n
		s1 = s1 + x(i)
		ss1 = ss1 + x(i)*x(i)
		s2 = s2 + y(i)
		ss2 = ss2 + y(i)*y(i)
		ss12 = ss12 + x(i)*y(i)
	ENDDO
	var1 = (ss1 - (s1*s1)/n) / (n-1)
	var2 = (ss2 - (s2*s2)/n) / (n-1)
	num = (ss12 - (s1*s2)/n) / (n-1)
	den = SQRT(var1*var2)

	IF( den > 0. )THEN
		f_val = num / den
	ELSE
		f_val = 0.
	ENDIF
END FUNCTION corr
