import casadi as cs
import torch
import l4casadi as l4c

if __name__ == '__main__':
    device = 'cpu'
    network = torch.load('./models/hcw_simple_model.pt').to(device)
    l4c_model = l4c.L4CasADi(network, has_batch=True, device='cpu')

    x_sym = cs.MX.sym('x', 9, 1)
    y_sym = l4c_model(x_sym)
    f = cs.Function('y', [x_sym], [y_sym])
    df = cs.Function('dy', [x_sym], [cs.jacobian(y_sym, x_sym)])
    # ddf = cs.Function('ddy', [x_sym], [cs.hessian(y_sym, x_sym)[0]])

    x = cs.DM([[0.], [2.], [3.], [1.], [2.], [3.], [0.], [0.01], [0.02]])
    print(l4c_model(x))
    print(f(x))
    print(df(x))
    # print(ddf(x))
