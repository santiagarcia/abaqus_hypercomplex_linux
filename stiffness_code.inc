      ! Generated Analytical Stiffness
      ! Common Subexpressions
      real(8) :: x0
      x0 = dN_u(1)**2
      real(8) :: x1
      x1 = dN_u(1)*u(1, 1) + dN_u(2)*u(7, 1) + 1
      real(8) :: x2
      x2 = x1**2
      real(8) :: x3
      x3 = dN_v(1)*u(2, 1) + dN_v(2)*u(6, 1) + dN_v(3)*u(8, 1) + dN_v(4)*u(12, 1)
      real(8) :: x4
      x4 = x3**2
      real(8) :: x5
      x5 = dN_v(1)*u(3, 1) - dN_v(2)*u(5, 1) + dN_v(3)*u(9, 1) - dN_v(4)*u(11, 1)
      real(8) :: x6
      x6 = x5**2
      real(8) :: x7
      x7 = x2 + x4 + x6
      real(8) :: x8
      x8 = sqrt(x7)
      real(8) :: x9
      x9 = 2*x8 - 2/x7
      real(8) :: x10
      x10 = 1d0/x8
      real(8) :: x11
      x11 = 1d0/(2*nu + 2)
      real(8) :: x12
      x12 = E_mod*w_det
      real(8) :: x13
      x13 = x11*x12
      real(8) :: x14
      x14 = (1.0d0/2.0d0)*Area*x10*x13
      real(8) :: x15
      x15 = x14*x9
      real(8) :: x16
      x16 = x7**(-3.0d0/2.0d0)
      real(8) :: x17
      x17 = (1.0d0/2.0d0)*Area*x13*x16*x9
      real(8) :: x18
      x18 = x17*x2
      real(8) :: x19
      x19 = dN_u(1)*x1
      real(8) :: x20
      x20 = 4/x7**2
      real(8) :: x21
      x21 = 2*x10
      real(8) :: x22
      x22 = x19*x20 + x19*x21
      real(8) :: x23
      x23 = x14*x19
      real(8) :: x24
      x24 = dN_v(1)*x3
      real(8) :: x25
      x25 = x17*x19
      real(8) :: x26
      x26 = x24*x25
      real(8) :: x27
      x27 = x20*x24 + x21*x24
      real(8) :: x28
      x28 = dN_v(1)*x5
      real(8) :: x29
      x29 = x25*x28
      real(8) :: x30
      x30 = x20*x28 + x21*x28
      real(8) :: x31
      x31 = dN_v(2)*x5
      real(8) :: x32
      x32 = x25*x31
      real(8) :: x33
      x33 = -x20*x31 - x21*x31
      real(8) :: x34
      x34 = dN_v(2)*x3
      real(8) :: x35
      x35 = x25*x34
      real(8) :: x36
      x36 = x20*x34 + x21*x34
      real(8) :: x37
      x37 = dN_u(2)*x1
      real(8) :: x38
      x38 = x20*x37 + x21*x37
      real(8) :: x39
      x39 = dN_u(1)*dN_u(2)
      real(8) :: x40
      x40 = x15*x39 - x18*x39
      real(8) :: x41
      x41 = dN_v(3)*x3
      real(8) :: x42
      x42 = x25*x41
      real(8) :: x43
      x43 = x20*x41 + x21*x41
      real(8) :: x44
      x44 = dN_v(3)*x5
      real(8) :: x45
      x45 = x25*x44
      real(8) :: x46
      x46 = x20*x44 + x21*x44
      real(8) :: x47
      x47 = dN_v(4)*x5
      real(8) :: x48
      x48 = x25*x47
      real(8) :: x49
      x49 = -x20*x47 - x21*x47
      real(8) :: x50
      x50 = dN_v(4)*x3
      real(8) :: x51
      x51 = x25*x50
      real(8) :: x52
      x52 = x20*x50 + x21*x50
      real(8) :: x53
      x53 = ddN_v(1)**2*x12
      real(8) :: x54
      x54 = dN_v(1)**2
      real(8) :: x55
      x55 = x15*x54
      real(8) :: x56
      x56 = x17*x54
      real(8) :: x57
      x57 = x14*x24
      real(8) :: x58
      x58 = x3*x5
      real(8) :: x59
      x59 = x56*x58
      real(8) :: x60
      x60 = x17*x24
      real(8) :: x61
      x61 = x31*x60
      real(8) :: x62
      x62 = ddN_v(1)*x12
      real(8) :: x63
      x63 = Iz*x62
      real(8) :: x64
      x64 = dN_v(1)*x15
      real(8) :: x65
      x65 = dN_v(2)*x64
      real(8) :: x66
      x66 = dN_v(1)*x17
      real(8) :: x67
      x67 = x4*x66
      real(8) :: x68
      x68 = -dN_v(2)*x67 + ddN_v(2)*x63 + x65
      real(8) :: x69
      x69 = x37*x60
      real(8) :: x70
      x70 = dN_v(3)*x64
      real(8) :: x71
      x71 = -dN_v(3)*x67 + ddN_v(3)*x63 + x70
      real(8) :: x72
      x72 = x44*x60
      real(8) :: x73
      x73 = x47*x60
      real(8) :: x74
      x74 = dN_v(4)*x64
      real(8) :: x75
      x75 = -dN_v(4)*x67 + ddN_v(4)*x63 + x74
      real(8) :: x76
      x76 = x14*x28
      real(8) :: x77
      x77 = Iy*x62
      real(8) :: x78
      x78 = -1.0d0/2.0d0*Area*E_mod*dN_v(1)*dN_v(2)*w_det*x11*x16*x6*x9 + ddN_v(2)* &
      x77 + x65
      real(8) :: x79
      x79 = x17*x37
      real(8) :: x80
      x80 = x28*x79
      real(8) :: x81
      x81 = -dN_v(3)*x6*x66 + ddN_v(3)*x77 + x70
      real(8) :: x82
      x82 = -1.0d0/2.0d0*Area*E_mod*dN_v(1)*dN_v(4)*w_det*x11*x16*x6*x9 + ddN_v(4)* &
      x77 + x74
      real(8) :: x83
      x83 = J_tor*x13
      real(8) :: x84
      x84 = x39*x83
      real(8) :: x85
      x85 = x14*x31
      real(8) :: x86
      x86 = ddN_v(2)**2*x12
      real(8) :: x87
      x87 = dN_v(2)**2
      real(8) :: x88
      x88 = x15*x87
      real(8) :: x89
      x89 = x17*x87
      real(8) :: x90
      x90 = x58*x89
      real(8) :: x91
      x91 = x31*x79
      real(8) :: x92
      x92 = x17*x31
      real(8) :: x93
      x93 = x41*x92
      real(8) :: x94
      x94 = ddN_v(2)*x12
      real(8) :: x95
      x95 = Iy*x94
      real(8) :: x96
      x96 = dN_v(2)*x15
      real(8) :: x97
      x97 = dN_v(3)*x96
      real(8) :: x98
      x98 = -1.0d0/2.0d0*Area*E_mod*dN_v(2)*dN_v(3)*w_det*x11*x16*x6*x9 + ddN_v(3)* &
      x95 + x97
      real(8) :: x99
      x99 = dN_v(4)*x96
      real(8) :: x100
      x100 = dN_v(2)*x17
      real(8) :: x101
      x101 = -dN_v(4)*x100*x6 + ddN_v(4)*x95 + x99
      real(8) :: x102
      x102 = x50*x92
      real(8) :: x103
      x103 = x14*x34
      real(8) :: x104
      x104 = x34*x79
      real(8) :: x105
      x105 = Iz*x94
      real(8) :: x106
      x106 = x100*x4
      real(8) :: x107
      x107 = -dN_v(3)*x106 + ddN_v(3)*x105 + x97
      real(8) :: x108
      x108 = -dN_v(4)*x106 + ddN_v(4)*x105 + x99
      real(8) :: x109
      x109 = x14*x37
      real(8) :: x110
      x110 = dN_u(2)**2
      real(8) :: x111
      x111 = x41*x79
      real(8) :: x112
      x112 = x44*x79
      real(8) :: x113
      x113 = x47*x79
      real(8) :: x114
      x114 = x50*x79
      real(8) :: x115
      x115 = x14*x41
      real(8) :: x116
      x116 = ddN_v(3)**2*x12
      real(8) :: x117
      x117 = dN_v(3)**2
      real(8) :: x118
      x118 = x117*x15
      real(8) :: x119
      x119 = x117*x17
      real(8) :: x120
      x120 = x119*x58
      real(8) :: x121
      x121 = x17*x41*x47
      real(8) :: x122
      x122 = ddN_v(3)*ddN_v(4)*x12
      real(8) :: x123
      x123 = dN_v(3)*dN_v(4)
      real(8) :: x124
      x124 = x123*x15
      real(8) :: x125
      x125 = Iz*x122 - x123*x17*x4 + x124
      real(8) :: x126
      x126 = x14*x44
      real(8) :: x127
      x127 = -1.0d0/2.0d0*Area*E_mod*dN_v(3)*dN_v(4)*w_det*x11*x16*x6*x9 + Iy*x122 + &
      x124
      real(8) :: x128
      x128 = x14*x47
      real(8) :: x129
      x129 = ddN_v(4)**2*x12
      real(8) :: x130
      x130 = dN_v(4)**2
      real(8) :: x131
      x131 = x130*x15
      real(8) :: x132
      x132 = x130*x17
      real(8) :: x133
      x133 = x132*x58
      real(8) :: x134
      x134 = x14*x50

      ! Stiffness Matrix Accumulation
      K_mat_local_val(1,1) = K_mat_local_val(1,1) + x0*x15 - x0*x18 + x22*x23
      K_mat_local_val(1,2) = K_mat_local_val(1,2) + (1.0d0/2.0d0)*Area*E_mod*dN_u(1)*w_det*x1*x10*x11*x27 - x26
      K_mat_local_val(1,3) = K_mat_local_val(1,3) + (1.0d0/2.0d0)*Area*E_mod*dN_u(1)*w_det*x1*x10*x11*x30 - x29
      K_mat_local_val(1,5) = K_mat_local_val(1,5) + x23*x33 + x32
      K_mat_local_val(1,6) = K_mat_local_val(1,6) + (1.0d0/2.0d0)*Area*E_mod*dN_u(1)*w_det*x1*x10*x11*x36 - x35
      K_mat_local_val(1,7) = K_mat_local_val(1,7) + x23*x38 + x40
      K_mat_local_val(1,8) = K_mat_local_val(1,8) + (1.0d0/2.0d0)*Area*E_mod*dN_u(1)*w_det*x1*x10*x11*x43 - x42
      K_mat_local_val(1,9) = K_mat_local_val(1,9) + (1.0d0/2.0d0)*Area*E_mod*dN_u(1)*w_det*x1*x10*x11*x46 - x45
      K_mat_local_val(1,11) = K_mat_local_val(1,11) + x23*x49 + x48
      K_mat_local_val(1,12) = K_mat_local_val(1,12) + (1.0d0/2.0d0)*Area*E_mod*dN_u(1)*w_det*x1*x10*x11*x52 - x51
