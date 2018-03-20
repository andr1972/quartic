//translated from //https://github.com/bmurray7/mersenne-twister-examples c source
using System;
using System.Threading.Tasks;

namespace Util
{
    class MT19937
    {
        // Assumes W = 32 (omitting this)
        private static ushort N = 624;
        private static uint M = 397;
        private static int R = 31;
        private static uint A = 0x9908B0DF;
        private static uint F = 1812433253;
        private static int U = 11;
        // Assumes D = 0xFFFFFFFF (omitting this)

        private static int S = 7;
        private static uint B = 0x9D2C5680;

        private static int T = 15;
        private static uint C = 0xEFC60000;

        private static int L = 18;

        private static uint MASK_LOWER = (uint)(1 << R) - 1;
        private static uint MASK_UPPER = (uint)(1 << R);

        private uint[] mt;
        private uint index;

        // Re-init with a given seed
        public MT19937(uint seed)
        {
            mt = new uint[N];
            mt[0] = seed;
            for (uint i = 1; i < N; i++)
            {
                mt[i] = (F * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i);
            }

            index = N;
        }

        void Twist()
        {
            for (uint i = 0; i < N; i++)
            {
                uint x = (mt[i] & MASK_UPPER) + (mt[(i + 1) % N] & MASK_LOWER);
                uint xA = x >> 1;
                if ((x & 0x1) != 0)
                    xA ^= A;
                mt[i] = mt[(i + M) % N] ^ xA;
            }
            index = 0;
        }

        // Obtain a 32-bit random number
        public uint Next()
        {
            uint i = index;

            if (index >= N)
            {
                Twist();
                i = index;
            }

            uint y = mt[i];
            index = i + 1;

            y ^= (y >> U);
            y ^= (y << S) & B;
            y ^= (y << T) & C;
            y ^= (y >> L);

            return y;
        }
    }
}
