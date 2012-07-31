    subroutine get_chi(w, t, f, ferr, n, order, B, chi, info)

        ! Get the ``chi`` vector for a dataset given a frequency
        ! and order.
        !
        ! :param w: The frequency of the lightcurve to fit.
        ! :param t: The vector of times in days.
        ! :param f: The vector of fluxes in arbitrary units.
        ! :param ferr: The vector of errors on the fluxes.
        ! :param n: The length of ``t``, ``f`` and ``ferr``.
        ! :param order: The integer order of model to fit.
        !
        ! :returns B: The best fit amplitude vector. This will
        !       have ``2 * order + 1`` values corresponding to the terms:
        !       ``(constant, sin(w t), cos(w t), sin(2 w t), ...)``.
        ! :returns chi: The residuals vector with ``n`` elements.
        ! :returns info: The result of the least squares fit: 0 for success.

        double precision, intent(in) :: w
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: t, f, ferr
        integer, intent(in) :: order

        double precision, dimension(1+2*order), intent(out) :: B
        double precision, dimension(n), intent(out) :: chi
        integer, intent(out) :: info

        double precision, dimension(n) :: ivar
        double precision, dimension(n, 1+2*order) :: A
        double precision, dimension(1+2*order, n) :: ATC
        double precision, dimension(1+2*order, 1+2*order) :: ATCA

        integer :: i, j, k, m, nlvl, liwork

        ! Constants.
        m = 2 * order + 1
        nlvl = max(0, int(log(real(m / 2))) + 1)

        ! Build the ``A`` matrix.
        A(:,1) = 1.0
        do i=1,order
            do j=1,n
                A(j,2*i)   = sin(i * t(j) * w)
                A(j,2*i+1) = cos(i * t(j) * w)
            enddo
        enddo

        ivar = 1.0 / ferr / ferr

        do i=1,m
            ATC(i,:) = A(:,i) * ivar
            ATCA(i,:) = 0.0
            do k=1,n
                ATCA(i,:) = ATCA(i,:) + ATC(i,k) * A(k,:)
            enddo
        enddo

        B(:) = 0.0
        do k=1,n
            B(:) = B(:) + ATC(:,k) * f(k)
        enddo

        lwork = -1
        liwork = 3 * m * nlvl + 11 * m

        call solve_square(m, ATCA, B, lwork, liwork, info)
        call solve_square(m, ATCA, B, lwork, liwork, info)

        if (info > 0) then
            write (*,*) "SVD didn't converge."
            return
        endif

        chi = f
        do i=1,m
            chi = chi - A(:,i) * B(i)
        enddo
        chi = chi * sqrt(ivar)

    end subroutine

    subroutine get_chi2(w, t, f, ferr, n, order, B, chi2, info)

        ! Get the sum of the squares of the residuals for a dataset given
        ! a frequency and order.
        !
        ! :param w: The frequency of the lightcurve to fit.
        ! :param t: The vector of times in days.
        ! :param f: The vector of fluxes in arbitrary units.
        ! :param ferr: The vector of errors on the fluxes.
        ! :param n: The length of ``t``, ``f`` and ``ferr``.
        ! :param order: The integer order of model to fit.
        !
        ! :returns B: The best fit amplitude vector. This will
        !       have ``2 * order + 1`` values corresponding to the terms:
        !       ``(constant, sin(w t), cos(w t), sin(2 w t), ...)``.
        ! :returns chi2: The chi-squared value for the dataset.
        ! :returns info: The result of the least squares fit: 0 for success.

        double precision, intent(in) :: w
        integer, intent(in) :: n
        double precision, dimension(n), intent(in) :: t, f, ferr
        integer, intent(in) :: order

        double precision, dimension(1+2*order), intent(out) :: B
        double precision, intent(out) :: chi2
        integer, intent(out) :: info

        integer :: k
        double precision, dimension(n) :: chi

        call get_chi(w, t, f, ferr, n, order, B, chi, info)

        chi2 = 0.0
        do k=1,n
            chi2 = chi2 + chi(k) * chi(k)
        enddo

    end subroutine

    subroutine solve_square(n, A, B, lwork, liwork, info)

        ! Lightweight wrapper around ``dgelsd`` from ``LAPACK`` specialized
        ! to square matrices.

        integer :: nrhs=1
        integer :: info, rank
        double precision :: rcond = -1.0

        integer :: n, lwork, liwork
        double precision, dimension(n,n) :: A
        double precision, dimension(n) :: B
        double precision, dimension(n) :: S

        double precision, dimension(max(1, lwork)) :: work
        integer, dimension(liwork) :: iwork

        S = 0.0
        work = 0.0
        iwork = 0
        info = 0

        call dgelsd(n, n, nrhs, A, n, B, n, S, rcond, rank, work, lwork, &
            iwork, info)

        lwork = int(work(1))

    end subroutine
