%%% @doc module enacl_ext implements various enacl extensions.
%%% <p>None of the extensions listed here are part of the official NaCl library.
%%% Functions may be removed without further notice if it suddenly ends up being
%%% better to do something differently than the solution given here.
%%% </p>
-module(enacl_ext).

-export([
         scramble_block_16/2,
         curve25519_scalarmult/2,
         curve25519_scalarmult_unrestrained/2,
         curve25519_recover_y/1,
         mod25519_sqrt/1
]).

%% @doc scramble_block_16/2 scrambles (encrypt) a block under a given key
%% The rules are that the block is 16 bytes and the key is 32 bytes. The block
%% is scrambled by means of the (secret) key. This makes it impossible for an
%% attacker to understand the original input for the scrambling. The intention
%% of this method is to protect counters from leaking to the outside world, by
%% scrambling them before they leave the system.
%%
%% Scrambling is done by means of the TEA algorithm (Tiny Encryption Algorithm)
%% It has known weaknesses and should probably not be used long-term going
%% forward, but CurveCP currently uses it for nonce scrambling.
%% @end
-spec scramble_block_16(binary(), binary()) -> binary().
scramble_block_16(Block, Key) ->
    enacl_nif:scramble_block_16(Block, Key).

%%% @doc
%%% Multiply a point by a scalar on curve25519, with Diffie-Helman
%%% appropriate masking of the scalar.
%%% @end
-spec curve25519_scalarmult(binary(), binary()) -> binary().
curve25519_scalarmult(N, Point) ->
    enacl_nif:curve25519_scalarmult(N, Point).

%%% @doc
%%% Given the X coordinate of a point on the curve, returns one of the two corresponding Y values (deterministically).
%%% Otherwise, returns {error, not_on_curve}.
%%% @end
curve25519_recover_y(N) ->
    enacl_nif:curve25519_ext_recover_y(N).

%%% @doc
%%% Given a number in the prime field (encoded as 256-bit little-endian binary), returns a square root (chosen deterministically).
%%% If no square root exists, returns {error, no_sqrt}.
%%% @end
mod25519_sqrt(N) ->
    enacl_nif:mod25519_sqrt(N).

%%% @doc
%%% Multiply a point by a scalar on curve25519, without masking of the scalar.
%%% @end
curve25519_scalarmult_unrestrained(N, Point) ->
    enacl_nif:curve25519_ext_scalarmult_unrestrained(N, Point).
