function f = objsppA0p(x)
if strcmp(x,'init')
    f.LB = [12.66 20.24 18.37 6.8 6.49 13.62 13.12 7.8 7.5 36.12 1.17 3.97 -77.99 -74.54 -70.67 -57.98 -68.22 -72.33 -66.77 -77.07 -73.86 -76.42 -74.37 -72.96 12.66 11.51 1.32 3.04 14.38 3.33 8.22 7.76 7.5 17.19 9.97 7.12 -79.04 -71.9 -66.69 -65.85 -75.68 -71.8 -68.52 -77.07 -74.75 -57.96 -72.61 -69.74 31.81 13.01 7.55 6.8 6.49 3.33 8.22 3.6 7.5 22.84 9.97 3.97 -77.99 -70.27 -70.67 -65.85 -75.68 -72.33 -65.84 -77.07 -56.56 -68.95 -74.37 -72.96 17 11.51 1.32 3.04 34.18 3.33 14.54 16.58 23.51 17.19 0.62 8.96 -79.04 -78.31 -66.69 -42.58 -77.08 -71.8 -68.52 -77.07 -74.75 -57.96 -72.61 -69.74 17 8.25 1.32 3.04 43.33 13.62 14.54 21.04 35.89 36.77 0.62 26.36 -79.04 -78.31 -58.93 -57.98 -77.08 -71.8 -44.87 -74.71 -74.75 -56.08 -56.08 -58.14 17 8.25 4.88 6.8 7.49 3.33 15.49 2.98 7.5 10.67 0.62 8.96 -77.82 -78.31 -66.69 -65.67 -77.08 -63.88 -69.83 -77.07 -73.86 -76.42 -76.47 -77.38 31.08 24.38 7.55 25.08 18.2 52.77 13.17 3.6 11.55 10.67 7.29 15.36 -59.74 -72.35 -57.77 -46.74 -52.58 -66.78 -25.05 -39.56 -65.7 -32.42 -50.79 -69.88 31.08 11.51 1.32 3.04 6.49 52.77 14.54 7.8 40.06 32.42 7.29 3.97 -79.04 -69.55 -70.67 -54.29 -66.51 -72.33 -63.73 -77.07 -74.75 -57.96 -74.37 -48.15 12.66 9.55 7.09 6.8 7.49 4.01 8.22 2.98 10.63 10.67 9.97 7.12 -77.82 -72.35 -66.69 -65.85 -75.68 -66.78 -69.83 -77.07 -65.7 -58.79 -76.47 -77.38 17 9.55 4.88 12.67 6.49 29.91 35.64 2.98 29.77 25.1 0.62 3.97 -77.99 -78.31 -70.67 -65.67 -77.08 -72.33 -69.83 -68.32 -73.86 -76.42 -76.47 -77.38 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    f.UB = [77.99 74.54 70.67 57.98 68.22 72.33 66.77 77.07 73.86 76.42 74.37 72.96 -12.66 -20.24 -18.37 -6.8 -6.49 -13.62 -13.12 -7.8 -7.5 -36.12 -1.17 -3.97 79.04 71.9 66.69 65.85 75.68 71.8 68.52 77.07 74.75 57.96 72.61 69.74 -12.66 -11.51 -1.32 -3.04 -14.38 -3.33 -8.22 -7.76 -7.5 -17.19 -9.97 -7.12 77.99 70.27 70.67 65.85 75.68 72.33 65.84 77.07 56.56 68.95 74.37 72.96 -31.81 -13.01 -7.55 -6.8 -6.49 -3.33 -8.22 -3.6 -7.5 -22.84 -9.97 -3.97 79.04 78.31 66.69 42.58 77.08 71.8 68.52 77.07 74.75 57.96 72.61 69.74 -17 -11.51 -1.32 -3.04 -34.18 -3.33 -14.54 -16.58 -23.51 -17.19 -0.62 -8.96 79.04 78.31 58.93 57.98 77.08 71.8 44.87 74.71 74.75 56.08 56.08 58.14 -17 -8.25 -1.32 -3.04 -43.33 -13.62 -14.54 -21.04 -35.89 -36.77 -0.62 -26.36 77.82 78.31 66.69 65.67 77.08 63.88 69.83 77.07 73.86 76.42 76.47 77.38 -17 -8.25 -4.88 -6.8 -7.49 -3.33 -15.49 -2.98 -7.5 -10.67 -0.62 -8.96 59.74 72.35 57.77 46.74 52.58 66.78 25.05 39.56 65.7 32.42 50.79 69.88 -31.08 -24.38 -7.55 -25.08 -18.2 -52.77 -13.17 -3.6 -11.55 -10.67 -7.29 -15.36 79.04 69.55 70.67 54.29 66.51 72.33 63.73 77.07 74.75 57.96 74.37 48.15 -31.08 -11.51 -1.32 -3.04 -6.49 -52.77 -14.54 -7.8 -40.06 -32.42 -7.29 -3.97 77.82 72.35 66.69 65.85 75.68 66.78 69.83 77.07 65.7 58.79 76.47 77.38 -12.66 -9.55 -7.09 -6.8 -7.49 -4.01 -8.22 -2.98 -10.63 -10.67 -9.97 -7.12 77.99 78.31 70.67 65.67 77.08 72.33 69.83 68.32 73.86 76.42 76.47 77.38 -17 -9.55 -4.88 -12.67 -6.49 -29.91 -35.64 -2.98 -29.77 -25.1 -0.62 -3.97 97 91 166 166 158 91 94 138 107 240 215 138 113 137 217 217 14 150 94 113 150 158 91 113 161 107 161 14 97 91 118 118 118 118 118 14 59 59 59 59 59 59 14 177 184 219 107 183 219 14 91 169 113 169 169 138 177 113 231 268 66 66 66 69 69 69 69 69 69 147 138 119 177 62 62 62 62 62 62 62 62 97 138 137 176 178 178 14 9 9 9 9 9 9 9 158 94 195 183 247 14 96 96 91 96 96 96 96 96 14 175 63 63 63 63 63 63 97 97 97 97 97 119 158 158 158 91 91 91 91 91 94 94 94 94 94 147 147 147 147 147 14 138 138 138 138 14 107 177 177 177 177 177 14 176 184 184 184 184 14 113 113 113 107 113 113 113 107 137 137 14];
    i = [1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 8 8 8 9 9 9 9 9 10 10 10 10 10 11 11 11 12 12 12 12 12 12 13 13 13 13 14 14 14 14 14 14 14 14 15 15 15 15 15 15 15 16 16 16 16 16 16 16 17 17 17 17 17 17 18 18 18 18 18 18 18 18 18 19 20 20 20 20 20 20 21 21 21 21 21 21 21 22 22 22 22 22 22 22 22 23 23 23 23 23 23 23 24 24 24 24 24 24 25 25 25 25 26 26 26 26 26 26 26 26 26 27 27 27 28 28 28 28 29 29 29 29 29 29 29 29 29 30 30 30 30 30 31 31 31 31 31 31 31 31 31 32 32 32 32 32 32 33 33 33 33 33 34 34 34 34 34 34 34 34 35 35 35 35 35 36 36 36 36 36 37 37 37 37 37 37 38 38 38 38 38 38 38 38 38 39 39 39 39 39 39 40 40 40 40 40 40 41 41 41 41 41 41 41 41 41 42 42 42 42 42 42 42 43 43 43 43 43 43 43 43 44 44 44 44 44 44 44 45 45 45 45 45 45 45 45 45 45 45 45 45];
    j = [241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 241 269 301 304 314 322 342 245 258 262 277 302 315 336 343 242 246 263 270 291 316 344 247 259 278 317 337 352 279 305 310 353 248 252 296 311 318 323 329 345 354 284 292 297 271 280 285 319 253 260 264 293 298 303 306 320 330 254 272 324 331 355 243 325 332 356 367 377 383 395 401 265 286 321 333 368 402 307 338 369 372 403 249 266 287 346 358 388 404 408 261 288 339 373 389 267 308 374 390 396 273 281 312 357 363 370 255 289 294 326 334 359 378 384 391 256 274 360 379 397 405 250 282 347 361 375 406 295 299 313 327 362 385 392 398 409 309 335 340 348 364 380 399 244 300 365 371 376 381 386 407 251 275 349 351 366 393 410 257 268 276 283 290 328 341 350 382 387 394 400 411];
    s = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    f.Aineq = full(sparse(i,j,s,45,411));
    f.bineq = [166 240 217 150 161 118 59 219 169 273 66 69 177 62 178 9 302 96 175 63 97 158 91 94 147 138 177 184 113 137 176 265 195 107 183 254 119 231 265 250 231 247 268 215 14]';
    i = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10];
    j = [241 269 301 304 314 322 342 358 359 360 361 362 245 258 262 277 302 315 336 343 363 364 365 366 242 246 263 270 291 316 344 367 368 369 370 371 247 259 278 317 337 352 372 373 374 375 376 279 305 310 353 377 378 379 380 381 382 248 252 296 311 318 323 329 345 354 383 384 385 386 387 284 292 297 388 389 390 391 392 393 394 271 280 285 319 395 396 397 398 399 400 253 260 264 293 298 303 306 320 330 401 402 403 404 405 406 407 254 272 324 331 355 408 409 410 411];
    s = [1 1 1 1 1 1 1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 -1 -1 -1 -1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1 -1 -1 -1 1 1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 1 1 1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 -1 -1 -1 -1];
    f.Aeq = full(sparse(i,j,s,10,411));
    f.beq = [0 0 0 0 0 0 0 0 0 0]';
    f.fitnessfcn = @objsppA0p;
    f.nonlcon = 'nlcsppA0p';
    f.nvars = 411;
    f.options.PopulationSize = 50;
    f.options.Generations = 10;
    f.options.ConstrBoundary = 'absorb';
else
    f = +21*x(241)+21*x(242)-27*x(243)-21*x(244)+43*x(245)+43*x(246)+43*x(247)+43*x(248)+3*x(249)-5*x(250)+2*x(251)+28*x(252)+28*x(253)+28*x(254)-21*x(255)-13*x(256)-12*x(257)+37*x(258)+37*x(259)+37*x(260)-13*x(261)+50*x(262)+50*x(263)+50*x(264)+6*x(265)+10*x(266)+2*x(267)+10*x(268)+13*x(269)+13*x(270)+13*x(271)+13*x(272)-37*x(273)-28*x(274)-28*x(275)-27*x(276)+14*x(277)+14*x(278)+14*x(279)+14*x(280)-36*x(281)-34*x(282)-26*x(283)+25*x(284)+25*x(285)-19*x(286)-15*x(287)-25*x(288)-24*x(289)-15*x(290)+31*x(291)+31*x(292)+31*x(293)-18*x(294)-10*x(295)+14*x(296)+14*x(297)+14*x(298)-27*x(299)-28*x(300)+36*x(301)+36*x(302)+36*x(303)+31*x(304)+31*x(305)+31*x(306)-19*x(307)-17*x(308)-9*x(309)+37*x(310)+37*x(311)-13*x(312)-4*x(313)+25*x(314)+25*x(315)+25*x(316)+25*x(317)+25*x(318)+25*x(319)+25*x(320)-19*x(321)+50*x(322)+50*x(323)+50*x(324)+2*x(325)+1*x(326)+9*x(327)+10*x(328)+43*x(329)+43*x(330)+43*x(331)-5*x(332)-1*x(333)-6*x(334)+3*x(335)+26*x(336)+26*x(337)-24*x(338)-24*x(339)-14*x(340)-14*x(341)+29*x(342)+29*x(343)+29*x(344)+29*x(345)-11*x(346)-19*x(347)-11*x(348)-12*x(349)-11*x(350)-7*x(351)+14*x(352)+14*x(353)+14*x(354)+14*x(355)-34*x(356)-36*x(357)-40*x(358)-49*x(359)-41*x(360)-48*x(361)-41*x(362)-50*x(363)-40*x(364)-42*x(365)-41*x(366)-48*x(367)-44*x(368)-50*x(369)-50*x(370)-42*x(371)-50*x(372)-50*x(373)-48*x(374)-48*x(375)-42*x(376)-48*x(377)-49*x(378)-41*x(379)-40*x(380)-42*x(381)-40*x(382)-48*x(383)-49*x(384)-41*x(385)-42*x(386)-40*x(387)-40*x(388)-50*x(389)-48*x(390)-49*x(391)-41*x(392)-41*x(393)-40*x(394)-48*x(395)-48*x(396)-41*x(397)-41*x(398)-40*x(399)-40*x(400)-48*x(401)-44*x(402)-50*x(403)-40*x(404)-41*x(405)-48*x(406)-42*x(407)-40*x(408)-41*x(409)-41*x(410)-40*x(411);
end