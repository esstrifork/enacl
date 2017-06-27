-module(add_or_subtract_test).

-include_lib("eunit/include/eunit.hrl").

%% Curve25519 definitions:
-define(PRIME25519, ((1 bsl 255) - 19) ).

identical_points_small_cases_test() ->
    %% Small X coordinates:
    lists:foreach(fun(I) ->
                          P = wrap(I),
                          R = enacl_ext:curve25519_add_or_subtract(P, P),
                          ?assertEqual({error, operator_failure}, R)
                  end,
                  lists:seq(0,100)).

identical_points_random_test() ->
    %% Random X coordinates:
    lists:foreach(fun(_) ->
                          P = wrap(rand_scalar()),
                          R = enacl_ext:curve25519_add_or_subtract(P, P),
                          ?assertEqual({error, operator_failure}, R)
                  end,
                  lists:seq(0,100)).

crypto_identical_points_test() ->
    %% Hidden identity:
    lists:foreach(fun(I) ->
                          X1 = I,
                          X2 = X1 + ?PRIME25519,
                          P = wrap(X1),
                          Q = wrap(X2),

                          R = enacl_ext:curve25519_add_or_subtract(P, Q),
                          ?assertEqual({error, operator_failure}, R)
                  end,
                  lists:seq(0, (1 bsl 255) - 1 - ?PRIME25519)).


small_cases_test() ->
    G = wrap(9),
    [check_scalar_case(G, NP, NQ)
     || NP <- lists:seq(0,15),
        NQ <- lists:seq(0,15)].

small_cases_other_bases_test() ->
    G = wrap(9),
    OtherBase = enacl_ext:curve25519_scalarmult_unrestrained(wrap(rand_scalar()), G),
    [check_scalar_case(OtherBase, NP, NQ)
     || NP <- lists:seq(0,15),
        NQ <- lists:seq(0,15)].

random_cases_test() ->
    [check_point_case(wrap(rand_scalar()), wrap(rand_scalar()))
     || _ <- lists:seq(1,200)].

special_cases_test() ->
    SpecialPoints = [0,1,2, 18,19,20] ++ lists:seq(?PRIME25519-2, (1 bsl 255)-1),
    [check_point_case(wrap(A), wrap(B))
     || A <- SpecialPoints,
        B <- SpecialPoints
    ].


check_scalar_case(Base, NP, NQ) ->
    P = enacl_ext:curve25519_scalarmult_unrestrained(wrap(NP), Base),
    Q = enacl_ext:curve25519_scalarmult_unrestrained(wrap(NQ), Base),

    R = enacl_ext:curve25519_add_or_subtract(P, Q),
    %% The operation should fail if-and-only-if either of:
    %% - P and Q are the same point
    %% - P is 0 or not on the curve
    %% - Q is 0 or not on the curve

    %% Cheap checks:
    ObviouslyInvalid =
        (unwrap(P) rem ?PRIME25519)==0 orelse
        (unwrap(Q) rem ?PRIME25519)==0 orelse
        ((unwrap(P)-unwrap(Q)) rem ?PRIME25519)==0,

    case ObviouslyInvalid of
        true ->
            ?assertEqual({error, operator_failure}, R);
        false ->
            case R of
                {error, operator_failure} ->
                    %% Less cheap checks:
                    YP = enacl_ext:curve25519_recover_y(P),
                    YQ = enacl_ext:curve25519_recover_y(Q),
                    ?assert(YP=={error, not_on_curve} orelse
                            YQ=={error, not_on_curve});
                _ when is_binary(R) ->
                    ?assertEqual(R, wrap(enacl_reference:add_or_subtract(unwrap(P), unwrap(Q)))),
                    ?assert(
                       R == enacl_ext:curve25519_scalarmult_unrestrained(wrap(NP + NQ), Base)
                       orelse
                       R == enacl_ext:curve25519_scalarmult_unrestrained(wrap(abs(NP - NQ)), Base)
                      )
            end
    end.

check_point_case(P, Q) ->
    R = enacl_ext:curve25519_add_or_subtract(P, Q),
    %% The operation should fail if-and-only-if either of:
    %% - P and Q are the same point
    %% - P is 0 or not on the curve
    %% - Q is 0 or not on the curve


    %% Cheap checks:
    ObviouslyInvalid =
        (unwrap(P) rem ?PRIME25519)==0 orelse
        (unwrap(Q) rem ?PRIME25519)==0 orelse
        ((unwrap(P)-unwrap(Q)) rem ?PRIME25519)==0,

    case ObviouslyInvalid of
        true ->
            ?assertEqual({error, operator_failure}, R);
        false ->
            case R of
                {error, operator_failure} ->
                    %% Less cheap checks:
                    YP = enacl_ext:curve25519_recover_y(P),
                    YQ = enacl_ext:curve25519_recover_y(Q),
                    ?assert(YP=={error, not_on_curve} orelse
                            YQ=={error, not_on_curve});
                _ when is_binary(R) ->
                    ?assertEqual(R, wrap(enacl_reference:add_or_subtract(unwrap(P), unwrap(Q))))
            end
    end.



wrap(X) -> <<X:256/little>>.
unwrap(<<X:256/little>>) -> X.

rand_scalar() ->
    <<X:256/little>> = crypto:rand_bytes(32),
    X rem ?PRIME25519.


timing_test() ->
    P1 = wrap(9),
    P2 = enacl_ext:curve25519_scalarmult_unrestrained(wrap(123), P1),
    report_time("scalarmult", fun() -> enacl_ext:curve25519_scalarmult(P1,P2) end,
                1000),
    report_time("scalarmult_unrestrained", fun() -> enacl_ext:curve25519_scalarmult_unrestrained(P1,P2) end,
                100),
    report_time("add_or_subtract", fun() -> enacl_ext:curve25519_add_or_subtract(P1,P2) end,
                1000),
    report_time("sqrt", fun() -> enacl_ext:mod25519_sqrt(P1) end,
                1000),
    report_time("recover_y", fun() -> enacl_ext:curve25519_recover_y(P2) end,
                1000),
    ok.

report_time(Msg, Fun, Iterations) when is_integer(Iterations) ->
    {T, _} = timer:tc(fun() -> timer_loop(Iterations, Fun) end),
    Seconds = (T*1.0e-6) / Iterations,
    Millis = Seconds*1.0e3,
    io:format(user, "Time for ~-30s: ~8.3f ms (average over ~b rounds)\n",
              [Msg, Millis, Iterations]).

timer_loop(0, _) -> ok;
timer_loop(I, Fun) ->
    Fun(),
    timer_loop(I-1, Fun).
