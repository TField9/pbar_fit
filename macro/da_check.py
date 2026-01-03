import inspect
from iminuit.cost import poisson_chi2
from iminuit.cost import template_chi2_da
print(inspect.getsource(template_chi2_da))
print("poisson_chi2函数的源代码定义：")
print(inspect.getsource(poisson_chi2))