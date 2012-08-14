subroutine dmim(sdat, nrs, ncs, ansm)
    use iface
    implicit none

    ! Input variables
    integer, intent(in) :: nrs, ncs
    integer, dimension(nrs, ncs), intent(in) :: sdat

    ! Output matrix
    real(kind=rdble), dimension(ncs, ncs), intent(out) :: ansm

    ! Local variables
    integer :: i, j, nok
    logical, dimension(nrs) :: ok
    ! Arrays to hold non-missing observations only
    integer, dimension(nrs) :: cvec, svec

    ! C function to get R integer NA value
    integer :: rnaint
    ! Local variable to hold R NA value
    integer :: naint

    ! Assign R NA value
    naint = rnaint()

    !$omp parallel do default(none) shared(ncs, sdat, naint, ansm) &
    !$omp private(i, j, ok, nok, cvec, svec) &
    !$omp schedule(dynamic)
    do i = 1, ncs
        do j = i, ncs
                ! Delete observations pairwise if either is missing
                ok = sdat(:,i) /= naint .and. sdat(:,j) /= naint

                nok = count(ok)

                cvec = pack(sdat(:,i), mask = ok)
                svec = pack(sdat(:,j), mask = ok)

                call dmi(cvec(1:nok), nok, svec(1:nok), nok, ansm(i,j))

                ! Symmetrise
                ansm(j,i) = ansm(i,j)
        end do
    end do
    !$omp end parallel do
end subroutine

subroutine dmi(v1, l1, v2, l2, ans) 
    use iface
    implicit none

    ! Input variables
    integer, intent(in) :: l1, l2
    integer, dimension(l1), intent(in) :: v1
    integer, dimension(l2), intent(in) :: v2

    ! Output variable
    real(kind=rdble), intent(out) :: ans

    ! Local variables
    integer i, j, lev1, lev2, k, l
    real(kind=rdble) :: tot

    ! Local dynamic arrays
    integer, dimension(:,:), allocatable :: tab
    real(kind=rdble), dimension(:,:), allocatable :: ptab
    real(kind=rdble), dimension(:), allocatable :: rv, cv

    lev1 = maxval(v1)
    lev2 = maxval(v2)

    ! Initialise dynamic arrays
    allocate(tab(lev1, lev2))
    tab = 0
    allocate(ptab(lev1, lev2))
    allocate(rv(lev1))
    allocate(cv(lev2))

    ! Cross tab for joint distribution
    do i = 1, l1
        tab(v1(i), v2(i)) = tab(v1(i), v2(i)) + 1
    end do

    ! Marginal counts (sum over joint)
    cv = sum(tab, dim = 1)
    rv = sum(tab, dim = 2)
    tot = sum(tab)

    ! Probability tables
    ptab = tab / tot
    rv = rv / tot
    cv = cv / tot

    ans = 0.0

    ! Calculate MI
    do i = 1, lev1
        do j = 1, lev2
            if (ptab(i, j) > 0) then
                ans = ans + ptab(i,j) * log(ptab(i,j) / (rv(i) * cv(j))) 
            end if
        end do
    end do

    deallocate(tab, ptab, rv, cv)
end subroutine
