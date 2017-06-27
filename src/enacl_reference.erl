-module(enacl_reference).

%%% ************************************************************
%%% * PURPOSE: Reference implementation of certain operations. *
%%% ************************************************************

-compile(export_all). % TODO

%% Curve25519 definitions:
-define(PRIME25519, ((1 bsl 255) - 19) ).
-define(CURVE_A, 486662).
-define(CURVE_Q, ((1 bsl 253) - 5) ).

dump(Msg,Num0) ->
    Num = normalize(Num0),
    true = (Num>=0),
    Bin = <<Num:256/little>>,
    io:format("ERLANG: ~s: ~p\n", [Msg, Bin]).

add_or_subtract(X1, X2) ->
    Y1 = recover_y_coordinate(X1),
    Y2 = recover_y_coordinate(X2),
    %% Formula: X3 = [(y1-y2)/(x1-x2)]² - A - x1 - x2.
    DX = X1-X2,
    DY = Y1-Y2,
    Slope = (DY * inv(DX, ?PRIME25519)) rem ?PRIME25519,
    SlopeSqr = (Slope * Slope) rem ?PRIME25519,
    X3 = (SlopeSqr - ?CURVE_A - (X1+X2)) rem ?PRIME25519,
    normalize(X3).

recover_y_coordinate(X) ->
    %% Curve equation: y² = x³ + Ax² + x.
    XSqr = X*X,
    YSqr = X*(XSqr + ?CURVE_A * X + 1),
    sqrt(YSqr).

normalize(X) ->
    (X + ?PRIME25519) rem ?PRIME25519.

%%% Calculates a square root in the prime field. Return 0 when there is no square root.
sqrt(X) ->
    %% Calculate by Tonelli-Shanks:
    SqrtExp = (?CURVE_Q+1) div 2,
    SqrtCand = powmod(X, SqrtExp, ?PRIME25519),
    Check = powmod(X, ?CURVE_Q, ?PRIME25519),
    InvDivisor = sqrt_divisor(Check),
    (InvDivisor * SqrtCand) rem ?PRIME25519.

%%% Table of (b² -> b^(-1)) for the possible results of (y²)^Q.
%%% The table contains: 0, 1, -1 (other sqrt of 1); no 4th roots of 1 exist.
sqrt_divisor(0) -> 0;
sqrt_divisor(1) -> 1;
sqrt_divisor(?PRIME25519-1) ->
    %% b²=-1; b found by Tonelli-Shanks algorithm.
    %% (b = 19681161376707505956807079304988542015446066515923890162744021073123829784752)
    38214883241950591754978413199355411911188925816896391856984770930832735035197;
sqrt_divisor(19681161376707505956807079304988542015446066515923890162744021073123829784752) ->
    0; % sqrt(-1)^(-1) => no square root.
sqrt_divisor(38214883241950591754978413199355411911188925816896391856984770930832735035197) ->
    0. % -sqrt(-1)^(-1) => no square root.



powmod(_A, 0, _N) -> 1;
powmod(A, 1, N) -> A rem N;
powmod(A, B, N) ->
    %% io:format("DB| powmod(~p,~p,~p)\n", [A, B, N]),
    A2 = powmod((A*A) rem N,
                B bsr 1, N),
    case B band 1 of
        0 ->
            A2;
        1 ->
            (A2*A) rem N
    end.

inv(X, P) ->
    powmod(X, P-2, P).

%%%========= Calculating square roots of elements of 2-power order:
tonelli_shanks(X) ->
    Z = 2,
    C = powmod(Z, ?CURVE_Q, ?PRIME25519),
    SqrtExp = (?CURVE_Q+1) div 2,
    M = 2, % Curve order = 2^M * Q

    R = powmod(X, SqrtExp, ?PRIME25519),
    T = powmod(X, ?CURVE_Q, ?PRIME25519),
    tonelli_shanks_loop(T, R, C, M).

tonelli_shanks_loop(0, R, _C, _M) -> R;
tonelli_shanks_loop(1, R, _C, _M) -> R;
tonelli_shanks_loop(T, R, C, M) ->
    % io:format("DB| tonelli_shanks_loop: T=~w R=~w C=~w M=~w\n", [T,R,C,M]),
    I = tonelli_shanks_count_squaring(T),
    (I<M) orelse error(no_sqrt),
    B = powmod(C, 1 bsl (M-I-1), ?PRIME25519),
    R2 = (R*B) rem ?PRIME25519,
    C2 = powmod(B, 2, ?PRIME25519),
    T2 = (T*C2) rem ?PRIME25519,
    tonelli_shanks_loop(T2, R2, C2, M).


tonelli_shanks_count_squaring(0) -> 0;
tonelli_shanks_count_squaring(1) -> 0;
tonelli_shanks_count_squaring(T) ->
    TSqr = (T*T) rem ?PRIME25519,
    1 + tonelli_shanks_count_squaring(TSqr).

%%%==================== Test: ========================================
-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").

-define(CURVE_ORDER,  ((1 bsl 252) + 16#14def9dea2f79cd65812631a5cf5d3ed) ).

sqrt_x2_returns_x_test() ->
    lists:foreach(fun(X) ->
                          X2 = sqrt(X*X),
                          ?assert(X2==X orelse X2==(?PRIME25519 - X))
                  end,
                  lists:seq(0,100)),
    lists:foreach(fun(_) ->
                          X = rand_scalar(),
                          X2 = sqrt(X*X),
                          ?assert(X2==X orelse X2==(?PRIME25519 - X))
                  end,
                  lists:seq(0,100)),
    ok.

recover_y_works_on_first_curve_points_test() ->
    lists:foreach(fun(I) ->
                          recover_y_works_on_curve_point(I)
                  end,
                  lists:seq(0,10)).

recover_y_works_on_random_curve_points_test() ->
    lists:foreach(fun(_) ->
                          recover_y_works_on_curve_point(rand_scalar())
                  end,
                  lists:seq(0,10)).

recover_y_works_on_curve_point(I) ->
    <<X:256/little>> = enacl_ext:curve25519_scalarmult_unrestrained(<<I:256/little>>, <<9:256/little>>),
    Y = recover_y_coordinate(X),
    Tmp = (X*X + ?CURVE_A*X + 1) rem ?PRIME25519,
    XX = X*Tmp rem ?PRIME25519,
    YY = (Y*Y) rem ?PRIME25519,
    ?assert(XX == YY).

sqrt_x2_does_not_break_test() ->
    lists:foreach(fun(X) ->
                          sqrt(X)
                  end,
                  lists:seq(0,500)).

sqrt_x2_does_not_break_random_test() ->
    lists:foreach(fun(_) ->
                          sqrt(rand_scalar())
                  end,
                  lists:seq(1,100)).

add_or_subtract_works_small_test() ->
    io:format("Making table...\n"),
    Table = lists:map(fun(I) ->
                      P_I = enacl_ext:curve25519_scalarmult_unrestrained(<<I:256/little>>, <<9:256/little>>),
                      <<X:256/little>> = P_I,
                      {I, X}
              end,
              lists:seq(0,10)),

    io:format("Creating test cases...\n"),
    TestPoints = [{M,N} || M <- lists:seq(0,10),
                           N <- lists:seq(0,10),
                           M =/= N,
                           M>0,
                           N>0,
                           M+N =< 10],

    io:format("Checking test cases...\n"),
    lists:foreach(fun({M,N}) ->
                          io:format("- Checking M=~p, N=~p...\n", [M,N]),
                          {_,XM} = lists:keyfind(M, 1, Table),
                          {_,XN} = lists:keyfind(N, 1, Table),
                          X = add_or_subtract(XM, XN),
                          {W,_} = lists:keyfind(X, 2, Table),
                          io:format("M=~p N=~p  => W=~p\n", [M,N,W]),
                          ?assert(W==M+N orelse W==abs(M-N))
                  end,
                  TestPoints).

add_or_subtract_works_random_test_() ->
    {timeout, 30, fun() ->
                          lists:foreach(fun(_) ->
                                                M = rand_scalar(),
                                                N = rand_scalar(),
                                                add_or_subtract_works(N,M)
                                        end,
                                        lists:seq(1,10))
                  end}.

add_or_subtract_works(N,M) ->
    io:format("Checking M=~p, N=~p...\n", [M,N]),
    CalcXForScalar = fun(I) ->
                             <<X:256/little>> = enacl_ext:curve25519_scalarmult_unrestrained(<<I:256/little>>, <<9:256/little>>),
                             io:format("CalcXForScalar: ~p -> X=~p\n", [I, X]),
                             X
                     end,
    XM = CalcXForScalar(M),
    XN = CalcXForScalar(N),
    X = add_or_subtract(XM, XN),
    io:format("- XM=~p  XN=~p  X=~p\n", [XM,XN,X]),

    ?assert(X == CalcXForScalar((M+N) rem ?CURVE_ORDER)
            orelse
            X == CalcXForScalar(abs(M-N))
            %% X == CalcXForScalar((M-N+?CURVE_ORDER) rem ?CURVE_ORDER)
            %% orelse
            %% X == CalcXForScalar((N-M+?CURVE_ORDER) rem ?CURVE_ORDER)
            ).


rand_scalar() ->
    <<X:256/little>> = crypto:rand_bytes(32),
    X rem ?PRIME25519.

-endif.
