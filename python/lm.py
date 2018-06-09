class LM():
    """
    Levenberg–Marquardt algorithm which allows you to specify which numpy package to use.
    For example, you can use the default numpy installation:

    import numpy
    lm = LM(numpy)

    or you can change the implementation to use the numpy version packaged with autograd:

    import autograd
    lmAG = LM(autograd.numpy)

    You might want to use this second version if you want, for instance, to be able to
    differentiate the output of LevenbergMarquardt, since autograd can automatically
    compute those derivatives for you.
    """

    def __init__(self, numpy):
        self.numpy = numpy

    def solve(self, f, df, x0, stop=None, step=1e-3, stepDown=None, stepUp=None, scaling=True):
        """
        This function uses the Levenberg–Marquardt algorithm to optimize a function f(x) given an initial guess x0.

        This is accomplished by repeatedly solving subproblems of the form
        (see p. 259 of "Numerical Optimization" by Nocedal and Wright):

                                                       2
             min  ||  [ df(x) / dx ]        [ f(x) ] ||
              dx  ||  [  step * I  ] * dx + [   0  ] ||

        to compute updates dx to the current best guess.
        If the new guess (x + dx) has a lower cost than the previous best guess, then the guess is updated.
        The main idea behind Levenberg–Marquardt is to adjust the step as we iterate,
        decreasing it when a new best guess is found
        (because this implies that we don't the regularization term),
        and increasing it when the new guess yields a higher cost than the initial guess
        (implying we need more regularization).
        The decrease and increase are controlled by the stepDown and stepUp coefficients, that is,
        step = step / stepDown or step = step * stepUp.

        TODO
        Note that the size of x and f(x) can have arbitrary dimensions. However, they must be the same type.
        If they weren't then we would have to do a cast when solving the subproblem. I believe that this
        would break the gradient propagation.
        """

        # Input handling
        if stepDown is None:
            stepDown = self.numpy.sqrt(10)
        if stepUp is None:
            stepUp = 3

        # If a stopping condition was not defined, define one
        if stop is None:
            maxIts = 20
            dcostTol = 1e-3
            dxTol = 1e-4

            def stop(x, fx, dfdx, cost, it, step, dcost, dx, scales, updated):
                if dcost is None or dx is None or self.numpy.isnan(dcost) or self.numpy.isinf(dcost):
                    return False
                if updated:
                    return it >= maxIts or self.numpy.abs(dcost) < dcostTol or self.numpy.linalg.norm(dx) < dxTol
                else:
                    return it >= maxIts

        # The subproblem requires a matrix of zeros that it doesn't make sense to keep recomputing
        xShape = self.numpy.shape(x0)
        dtype = x0.dtype
        xSize = self.numpy.size(x0)
        zeros = self.numpy.zeros((xSize, 1), dtype=dtype)

        def F(x):
            """ This function simply evaluates f(x) and then vectorizes it so that the output is (L,1). """
            return self.numpy.reshape(f(x), (-1, 1))

        defaultScales = self.numpy.ones((xSize,))

        def dF(x):
            """ The output fx of F(x) is always a matrix with a shape of the form (L,1).
            We need df to return something of size (L,n), where n is the size of x.
            Hence this function needs to compute the gradient, and then vectorize it. """
            dfdx = self.numpy.reshape(df(x), (-1, self.numpy.size(x)))
            if scaling:
                scales = self.numpy.maximum(1e-10, self.numpy.linalg.norm(dfdx, axis=0))
                dfdx /= scales
            else:
                scales = defaultScales
            return dfdx, scales

        def subproblem(dfdx, fx, step, scales):
            Istep = step * self.numpy.identity(xSize)
            A = self.numpy.vstack((dfdx, Istep))
            B = self.numpy.vstack((fx, zeros))
            dx = self.numpy.linalg.lstsq(A, B, rcond=None)[0]

            # We had to vectorize the problem so that the dimensions work for Levenburg Marquardt.
            # So we unvectorize here
            dx = dx.ravel() / scales
            dx = self.numpy.reshape(dx, xShape)
            return dx

        def loopBody(x0, fx0, dfdx0, cost0, it0, step0, dcost0, dx0, scales0, updated0):

            """
            Each time the body is evaluated, we attempt to step in a certain direction.
            If the step decreased the cost, then we accept the solution and update the step.
            Otherwise, we just modify the step.
            """
            # Propose an update based on the step
            dx = subproblem(dfdx0, fx0, step0, scales0)
            x = x0 - dx

            # Then compute its cost and jacobian
            fx = F(x)
            cost = self.numpy.linalg.norm(fx)

            # If the cost is lower, use the update
            dcost = cost0 - cost
            if self.numpy.isnan(dcost) or self.numpy.isnan(dcost):
                updated = False
            elif dcost > 0:
                updated = True
            else:
                updated = False

            print(it0, step0, cost0, cost, dcost, updated)

            if updated:
                dfdx, scales = dF(x)
                return x, fx, dfdx, cost, it0 + 1, step0 / stepDown, dcost, dx, scales, updated
            else:
                return x0, fx0, dfdx0, cost0, it0 + 1, step0 * stepUp, dcost0, dx0, scales0, updated

        # Compute the function, jacobian, and cost at the initial guess
        fx0 = F(x0)
        dfdx0, scales0 = dF(x0)
        cost0 = self.numpy.linalg.norm(fx0)
        it0 = 0
        step0 = step
        updated0 = False

        # Note that these values won't be used in the first iteration.
        # The important thing is that they are the right size and type
        dcost0 = None
        dx0 = None

        # Now, let's loop
        variables = x0, fx0, dfdx0, cost0, it0, step0, dcost0, dx0, scales0, updated0
        while not stop(*variables):
            variables = loopBody(*variables)
        return variables


LevenbergMarquardt = LM
