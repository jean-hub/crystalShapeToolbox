function hCost_13 = hessD2Mod_13_3(in1,in2,in3)
%HESSD2MOD_13_3
%    HCOST_13 = HESSD2MOD_13_3(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 5.9.
%    23-Sep-2014 14:06:30

d_131 = in3(1,:);
d_132 = in3(2,:);
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
yy_131 = in2(1,:);
yy_132 = in2(2,:);
yy_133 = in2(3,:);
t2 = q1.^2;
t3 = q2.^2;
t4 = q3.^2;
t5 = q4.^2;
t6 = q1.*q2.*2.0;
t39 = q1.*yy_132.*2.0;
t40 = q2.*yy_133.*2.0;
t41 = q4.*yy_131.*2.0;
t7 = t39-t40+t41;
t14 = q1.*1.224646799147353e-16;
t15 = q3.*2.0;
t16 = t14+t15;
t17 = t16.*yy_131;
t18 = q1.*2.0;
t19 = q3.*1.224646799147353e-16;
t20 = t18-t19;
t21 = t20.*yy_133;
t22 = q2.*2.0;
t23 = q4.*1.224646799147353e-16;
t24 = t22+t23;
t25 = t24.*yy_132;
t8 = -t17+t21+t25;
t9 = t2-t3+t4-t5;
t10 = q1.*q4.*2.0;
t11 = q2.*q3.*2.0;
t12 = t10+t11;
t13 = q3.*q4.*2.0;
t26 = q1.*q3.*1.224646799147353e-16;
t27 = q2.*q4.*1.224646799147353e-16;
t28 = -t2+t3+t4-t5+t26+t27;
t29 = q1.*q3.*2.0;
t30 = t2.*6.123233995736766e-17;
t31 = t3.*6.123233995736766e-17;
t32 = t4.*6.123233995736766e-17;
t33 = t5.*6.123233995736766e-17;
t46 = q2.*q4.*2.0;
t34 = t29+t30+t31-t32-t33-t46;
t35 = q1.*q4.*1.224646799147353e-16;
t48 = q2.*q3.*1.224646799147353e-16;
t36 = t6+t13+t35-t48;
t37 = t36.*yy_132;
t45 = t28.*yy_133;
t47 = t34.*yy_131;
t38 = d_131+t37-t45-t47;
t42 = q2.*1.224646799147353e-16;
t43 = q4.*2.0;
t44 = t42-t43;
t49 = t6-t13;
t50 = t49.*yy_133;
t52 = t9.*yy_132;
t53 = t12.*yy_131;
t51 = d_132+t50-t52-t53;
t54 = t51.*yy_133.*(1.0./2.0);
t55 = t44.*yy_131;
t56 = t24.*yy_133;
t65 = t20.*yy_132;
t57 = t55+t56-t65;
t58 = t38.*yy_132.*(1.0./2.0);
t59 = q1.*yy_133.*2.0;
t60 = q2.*yy_132.*2.0;
t64 = q3.*yy_131.*2.0;
t61 = t59+t60-t64;
t62 = t54+t58-t8.*t57.*(1.0./4.0)-t7.*t61.*(1.0./4.0);
t63 = yy_131.*1.224646799147353e-16;
t66 = t20.*yy_131;
t67 = t16.*yy_133;
t68 = t44.*yy_132;
t69 = t66+t67+t68;
t70 = t38.*yy_132.*3.061616997868383e-17;
t71 = q2.*yy_131.*2.0;
t72 = q3.*yy_132.*2.0;
t73 = q4.*yy_133.*2.0;
t74 = t71+t72+t73;
t75 = t16.*yy_132;
t76 = t24.*yy_131;
t93 = t44.*yy_133;
t77 = t75+t76-t93;
t78 = yy_131.*2.0;
t79 = yy_133.*1.224646799147353e-16;
t80 = q1.*yy_131.*2.0;
t81 = q3.*yy_133.*2.0;
t94 = q4.*yy_132.*2.0;
t82 = t80+t81-t94;
t83 = q2.*t51.*(1.0./2.0);
t84 = t20.*t38.*(1.0./4.0);
t85 = t24.*t38.*(1.0./4.0);
t86 = t78+t79;
t87 = t7.*t74.*(1.0./4.0);
t88 = t87-t8.*t69.*(1.0./4.0)-t38.*t86.*(1.0./4.0);
t89 = t57.*t69.*(1.0./4.0);
t95 = t51.*yy_131.*(1.0./2.0);
t90 = -t70+t89-t95-t61.*t74.*(1.0./4.0);
t91 = t51.*yy_132.*(1.0./2.0);
t92 = yy_133.*2.0;
t96 = t8.*t77.*(1.0./4.0);
t97 = t7.*t82.*(1.0./4.0);
t98 = t78-t79;
t99 = t38.*t98.*(1.0./4.0);
t100 = t99-t57.*t77.*(1.0./4.0)-t61.*t82.*(1.0./4.0);
t101 = t74.*t82.*(1.0./4.0);
t102 = -t54+t58+t101-t69.*t77.*(1.0./4.0);
t103 = t63+t92;
t104 = q1.*t51.*(1.0./2.0);
t105 = q4.*t51.*(1.0./2.0);
t106 = t16.*t38.*(1.0./4.0);
t107 = t7.*t12.*(1.0./4.0);
t108 = t34.*t57.*(1.0./4.0);
t118 = q3.*t51.*(1.0./2.0);
t119 = t38.*t44.*(1.0./4.0);
t109 = t108-t118-t119-t12.*t61.*(1.0./4.0);
t110 = t34.*t69.*(1.0./4.0);
t111 = t12.*t74.*(1.0./4.0);
t112 = -t83-t84+t110+t111;
t113 = t12.*t82.*(1.0./4.0);
t114 = t85-t104+t113-t34.*t77.*(1.0./4.0);
t115 = t8.*t36.*(1.0./4.0);
t116 = t7.*t9.*(1.0./4.0);
t117 = t83+t84-t9.*t61.*(1.0./4.0)-t36.*t57.*(1.0./4.0);
t120 = t9.*t74.*(1.0./4.0);
t121 = t36.*t77.*(1.0./4.0);
t122 = t9.*t82.*(1.0./4.0);
t123 = t105+t106+t121+t122;
t124 = t9.*t12.*(1.0./4.0);
t125 = t124-t34.*t36.*(1.0./4.0);
t126 = t83+t84-t8.*t28.*(1.0./4.0)-t7.*t49.*(1.0./4.0);
t127 = t28.*t57.*(1.0./4.0);
t128 = t49.*t61.*(1.0./4.0);
t129 = -t85+t104+t127+t128;
t130 = t28.*t69.*(1.0./4.0);
t131 = t28.*t34.*(1.0./4.0);
t132 = t131-t12.*t49.*(1.0./4.0);
t133 = t9.*t49.*(-1.0./4.0)-t28.*t36.*(1.0./4.0);
hCost_13 = reshape([yy_132.*(d_132-t9.*yy_132-t12.*yy_131+yy_133.*(t6-q3.*q4.*2.0)).*(-1.0./2.0)-t38.*(t63-yy_133.*2.0).*(1.0./4.0)+t7.^2.*(1.0./4.0)+t8.^2.*(1.0./4.0),t62,t88,t70-t95+t96+t97,-t105-t106+t107-t8.*t34.*(1.0./4.0),t85-t104+t115+t116,t126,t62,t91-t38.*t103.*(1.0./4.0)+t57.^2.*(1.0./4.0)+t61.^2.*(1.0./4.0),t90,t100,t109,t117,t129,t88,t90,-t91+t38.*(t63-t92).*(1.0./4.0)+t69.^2.*(1.0./4.0)+t74.^2.*(1.0./4.0),t102,t112,-t118-t119+t120-t36.*t69.*(1.0./4.0),-t105-t106+t130-t49.*t74.*(1.0./4.0),t70+t96+t97-t51.*yy_131.*(1.0./2.0),t100,t102,t91+t38.*t103.*(1.0./4.0)+t77.^2.*(1.0./4.0)+t82.^2.*(1.0./4.0),t114,t123,-t118-t119-t28.*t77.*(1.0./4.0)-t49.*t82.*(1.0./4.0),t107-q4.*t51.*(1.0./2.0)-t8.*t34.*(1.0./4.0)-t16.*t38.*(1.0./4.0),t109,t112,t114,t12.^2.*(1.0./4.0)+t34.^2.*(1.0./4.0),t125,t132,t85+t115+t116-q1.*t51.*(1.0./2.0),t117,t120-q3.*t51.*(1.0./2.0)-t38.*t44.*(1.0./4.0)-t36.*t69.*(1.0./4.0),t123,t125,t9.^2.*(1.0./4.0)+t36.^2.*(1.0./4.0),t133,t126,t129,t130-q4.*t51.*(1.0./2.0)-t16.*t38.*(1.0./4.0)-t49.*t74.*(1.0./4.0),q3.*t51.*(-1.0./2.0)-t38.*t44.*(1.0./4.0)-t28.*t77.*(1.0./4.0)-t49.*t82.*(1.0./4.0),t132,t133,t28.^2.*(1.0./4.0)+t49.^2.*(1.0./4.0)],[7, 7]);