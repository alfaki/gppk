function f = objsppB0(x)
if strcmp(x,'init')
    f.LB = zeros(1,558);
    f.UB = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 149 206 137 58 58 58 58 58 58 58 58 20 58 21 191 182 60 175 207 149 149 82 60 175 90 124 99 21 78 97 149 60 145 175 65 124 159 191 191 78 151 151 145 151 65 90 151 77 77 77 77 77 65 77 77 21 16 77 58 58 58 58 58 58 58 7 35 35 35 35 7 175 65 99 172 16 181 137 97 102 60 102 102 90 102 22 102 149 82 175 99 159 278 16 143 143 99 143 143 143 16 78 143 97 82 124 124 124 124 90 124 124 124 22 109 60 109 109 63 109 109 7 67 67 67 60 65 67 67 67 67 23 23 23 23 23 16 23 7 120 120 120 120 120 63 20 120 149 180 99 82 145 100 124 182 20 22 16 149 82 182 175 100 289 182 22 126 126 126 99 126 126 16 126 97 151 124 99 151 151 191 149 149 182 145 65 88 242 191 149 149 145 65 124 99 159 88 29 29 29 29 29 20 22 16 29 29 47 47 47 47 47 47 47 16 47 97 149 182 65 100 119 99 254 20 108 242 149 82 60 100 90 20 236 97 149 60 119 90 21 97 116 100 116 116 116 175 172 175 97 149 175 65 100 90 124 99 149 149 119 90 20 158 97 121 121 121 121 78 172 191 191 21 136 137 63 97 97 97 20 97 16 78 149 21 149 149 16 149 7 149 149 149 149 137 82 82 78 82 150 21 108 182 7 88 20 180 21 180 108 180 78 180 60 60 60 20 60 7 145 63 108 145 22 145 175 20 175 108 175 175 137 65 65 63 65 65 16 65 100 100 21 100 100 100 22 100 100 100 7 119 63 119 119 119 88 21 119 7 90 90 90 90 90 7 124 124 21 78 124 99 99 21 99 99 22 99 159 159 149 88 108 22 16];
    i = [1 1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 9 9 9 9 9 10 10 10 10 10 10 10 11 11 11 11 11 11 11 11 11 12 12 12 12 12 12 12 13 13 13 13 13 13 13 13 13 14 14 14 14 14 14 14 14 14 14 14 15 15 15 15 15 15 15 15 16 16 16 16 16 16 16 16 16 17 17 17 17 17 17 17 17 18 18 18 18 18 18 18 18 19 19 19 20 20 20 20 20 20 20 20 21 21 21 21 21 21 21 21 22 22 22 22 22 22 22 22 23 23 23 23 23 23 24 24 24 24 24 24 24 24 25 25 25 25 25 25 25 25 25 26 26 26 26 26 26 26 26 26 26 27 27 27 27 27 27 27 27 27 28 28 28 28 28 28 28 28 28 28 28 29 29 29 29 29 29 29 30 30 30 30 30 30 31 31 31 31 31 31 32 32 32 33 33 33 33 33 33 33 33 34 34 34 34 34 34 35 35 35 35 35 35 36 36 36 36 36 36 37 37 37 37 37 37 37 37 37 37 38 38 38 38 38 38 38 38 38 38 38 38 38 38 39 39 39 39 39 39 39 39 39 39 39 39 39 39 40 40 40 40 40 40 40 40 40 40 40 41 41 41 41 41 41 41 41 41 41 42 42 42 42 42 42 42 42 42 42 43 43 43 43 43 43 43 43 44 44 44 44 44 44 44 44 44 44 44 44 44 45 45 45 45 45 45 45 45 45 45 45 45 45 45 46 46 46 46 46 46 46 46 46 46 46 47 47 47 47 47 47 48 48 48 48 48 48 48 48 49 49 49 49 49 49 49 49 49 49 49 50 50 50 50 50 50 50 50 51 51 51 51 51 51 51 51 51 51 51 52 52 52 52 52 52 52 52 52 53 53 53 53 53 53 53 53 53 53 53 54 54 54 54 54 54 54 54 54 54 55 55 55 55 55 55 55 56 56 56 56 56 56 56 56 56 56 56 57 57 57 57 57 57 57 57 57 57 58 58 58 58 58 58 58 58 58 59 59 59 59 59 59 59 60 60 60 60 60 60 60 60 60 60 60 61 61 61 61 61 61 61 61 61 62 62 62 62 62 62 62 62 63 63 63 63 63 63 63 63 63 63 63 63 64 64 64 64 64 64 64 64 64 65 65 65 65 65 65 65 65 65 66 66 66 66 66 66 66 66 66 66 66 67 67 67 67 67 67 67 67 67 68 68 68 68 68 68 68 68 68 68 68 68 68 69 69 69 69 69 69 69 69 69 69 69 70 70 70 70 70 70 70 70 70 70 70 70 70 70 71 71 71 71 71 71 71 71 71 72 72 72 72 72 72 72 72 72 72 73 73 73 73 73 73 73 73 73 73];
    j = [175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 394 395 396 397 398 399 400 401 402 403 404 405 406 407 408 409 410 411 412 413 414 415 416 417 418 419 420 421 422 423 424 425 426 427 428 429 430 431 432 433 434 435 436 437 438 439 440 441 442 443 189 215 341 355 363 421 178 204 253 278 349 391 409 415 424 438 194 222 254 314 333 342 356 364 372 392 402 410 425 432 175 179 195 205 223 233 262 269 297 322 357 365 382 433 180 196 224 234 241 263 279 298 325 334 403 190 225 235 280 306 315 335 350 358 393 181 216 242 281 289 299 307 323 343 439 191 197 206 255 290 300 404 411 207 217 236 256 270 282 291 316 326 359 366 383 416 192 198 208 218 226 246 257 264 283 317 336 373 426 440 209 219 227 247 301 360 367 374 384 394 427 327 337 395 405 417 428 182 243 308 385 396 412 418 434 199 220 228 258 284 302 309 406 413 429 435 200 210 292 328 351 368 430 441 183 201 248 265 271 324 344 352 369 397 431 211 221 259 266 272 285 310 370 386 184 237 249 318 422 444 465 507 514 540 552 185 229 267 338 375 398 494 508 525 541 238 293 319 450 495 509 526 273 286 303 329 339 419 451 458 466 470 527 274 294 445 452 488 500 510 528 534 553 304 345 376 442 453 467 489 529 554 305 361 371 479 490 530 555 186 320 330 377 399 407 436 454 480 491 501 239 346 353 423 446 481 502 511 545 187 287 387 455 468 474 515 546 188 202 230 414 447 459 475 482 516 531 542 547 193 212 354 388 408 460 471 483 517 400 476 484 496 503 518 535 548 556 213 295 362 401 437 461 485 497 504 519 549 260 288 331 340 378 498 520 550 557 231 250 268 275 311 332 347 379 389 456 462 512 558 203 214 232 276 390 443 457 472 486 492 543 176 251 277 321 380 420 463 477 487 499 505 513 521 536 261 312 381 448 473 522 532 537 544 177 244 252 348 449 469 506 523 538 551 240 245 296 313 464 478 493 524 533 539];
    s = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    f.Aineq = full(sparse(i,j,s,73,558));
    f.bineq = [284 58 207 217 191 151 77 58 35 181 102 278 143 124 109 67 23 120 303 248 289 126 151 258 301 29 47 254 236 225 116 175 208 158 121 191 97 149 149 82 182 180 60 145 175 65 100 119 90 124 99 159 172 290 63 182 280 149 88 20 287 150 21 240 108 242 22 16 78 206 136 137 7]';
    i = [1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7 8 8 8 8 8 8 8 8 9 9 9 9 9 9 9 9 9 9 9 9 9 10 10 10 10 10 10 10 10 10 10 10 10 10 10 11 11 11 11 11 11 11 11 11 11 11 12 12 12 12 12 12 13 13 13 13 13 13 13 13 14 14 14 14 14 14 14 14 14 14 14 15 15 15 15 15 15 15 15 16 16 16 16 16 16 16 16 16 16 16 17 17 17 17 17 17 17 17 17];
    j = [8 28 107 115 121 158 2 20 52 68 111 137 149 154 159 171 12 35 53 91 102 108 116 122 129 138 144 150 160 167 1 3 13 21 36 42 59 64 80 95 117 123 132 168 4 14 37 43 46 60 69 81 98 103 145 9 38 44 70 86 92 104 112 118 139 5 29 47 71 76 82 87 96 109 172 10 15 22 54 77 83 146 151 23 30 45 55 65 72 78 93 99 119 124 133 155 11 16 24 31 39 49 56 61 73 94 105 130 161 173 25 32 40 50 84 120 125 131 134 140 162 100 106 141 147 156 163 6 48 88 135 142 152 157 169 17 33 41 57 74 85 89 148 153 164 170 18 26 79 101 113 126 165 174 7 19 51 62 66 97 110 114 127 143 166 27 34 58 63 67 75 90 128 136];
    s = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    f.Aeq = full(sparse(i,j,s,17,558));
    f.beq = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
    f.fitnessfcn = @objsppB0;
    f.nonlcon = 'nlcsppB0';
    f.nvars = 558;
    f.options.PopulationSize = 50;
    f.options.Generations = 10;
    f.options.ConstrBoundary = 'absorb';
else
    f = +45*x(175)+3*x(176)+2*x(177)+24*x(178)+24*x(179)+24*x(180)+24*x(181)+24*x(182)+24*x(183)-17*x(184)-26*x(185)-23*x(186)-17*x(187)-22*x(188)+17*x(189)+17*x(190)+17*x(191)+17*x(192)-32*x(193)+28*x(194)+28*x(195)+28*x(196)+28*x(197)+28*x(198)+28*x(199)+28*x(200)+28*x(201)-18*x(202)-20*x(203)+37*x(204)+37*x(205)+37*x(206)+37*x(207)+37*x(208)+37*x(209)+37*x(210)+37*x(211)-12*x(212)-9*x(213)-11*x(214)+45*x(215)+45*x(216)+45*x(217)+45*x(218)+45*x(219)+45*x(220)+45*x(221)+41*x(222)+41*x(223)+41*x(224)+41*x(225)+41*x(226)+41*x(227)+41*x(228)-9*x(229)-5*x(230)+1*x(231)-7*x(232)+33*x(233)+33*x(234)+33*x(235)+33*x(236)-8*x(237)-16*x(238)-10*x(239)-9*x(240)+46*x(241)+46*x(242)+46*x(243)+3*x(244)+4*x(245)+35*x(246)+35*x(247)+35*x(248)-6*x(249)-5*x(250)-7*x(251)-8*x(252)+47*x(253)+47*x(254)+47*x(255)+47*x(256)+47*x(257)+47*x(258)+47*x(259)-2*x(260)-1*x(261)+25*x(262)+25*x(263)+25*x(264)+25*x(265)+25*x(266)-25*x(267)-15*x(268)+27*x(269)+27*x(270)+27*x(271)+27*x(272)-14*x(273)-16*x(274)-13*x(275)-21*x(276)-15*x(277)+23*x(278)+23*x(279)+23*x(280)+23*x(281)+23*x(282)+23*x(283)+23*x(284)+23*x(285)-18*x(286)-18*x(287)-26*x(288)+20*x(289)+20*x(290)+20*x(291)+20*x(292)-29*x(293)-23*x(294)-26*x(295)-22*x(296)+23*x(297)+23*x(298)+23*x(299)+23*x(300)+23*x(301)+23*x(302)-18*x(303)-21*x(304)-20*x(305)+12*x(306)+12*x(307)+12*x(308)+12*x(309)+12*x(310)-28*x(311)-36*x(312)-30*x(313)+43*x(314)+43*x(315)+43*x(316)+43*x(317)+2*x(318)-6*x(319)-4*x(320)+1*x(321)+28*x(322)+28*x(323)+28*x(324)+42*x(325)+42*x(326)+42*x(327)+42*x(328)+1*x(329)-5*x(330)-7*x(331)+2*x(332)+27*x(333)+27*x(334)+27*x(335)+27*x(336)+27*x(337)-23*x(338)-14*x(339)-22*x(340)+41*x(341)+41*x(342)+41*x(343)+41*x(344)-3*x(345)-2*x(346)+1*x(347)-2*x(348)+21*x(349)+21*x(350)+21*x(351)+21*x(352)-22*x(353)-28*x(354)+19*x(355)+19*x(356)+19*x(357)+19*x(358)+19*x(359)+19*x(360)-24*x(361)-27*x(362)+26*x(363)+26*x(364)+26*x(365)+26*x(366)+26*x(367)+26*x(368)+26*x(369)+26*x(370)-17*x(371)+46*x(372)+46*x(373)+46*x(374)-4*x(375)+2*x(376)-1*x(377)-3*x(378)+6*x(379)+4*x(380)-2*x(381)+11*x(382)+11*x(383)+11*x(384)+11*x(385)+11*x(386)-30*x(387)-38*x(388)-29*x(389)-37*x(390)+43*x(391)+43*x(392)+43*x(393)+43*x(394)+43*x(395)+43*x(396)+43*x(397)-7*x(398)-4*x(399)+2*x(400)-3*x(401)+41*x(402)+41*x(403)+41*x(404)+41*x(405)+41*x(406)-6*x(407)-8*x(408)+28*x(409)+28*x(410)+28*x(411)+28*x(412)+28*x(413)-18*x(414)+39*x(415)+39*x(416)+39*x(417)+39*x(418)-2*x(419)-3*x(420)+16*x(421)-25*x(422)-27*x(423)+32*x(424)+32*x(425)+32*x(426)+32*x(427)+32*x(428)+32*x(429)+32*x(430)+32*x(431)+50*x(432)+50*x(433)+50*x(434)+50*x(435)+3*x(436)+4*x(437)+22*x(438)+22*x(439)+22*x(440)+22*x(441)-22*x(442)-26*x(443)-41*x(444)-43*x(445)-43*x(446)-46*x(447)-48*x(448)-43*x(449)-49*x(450)-41*x(451)-43*x(452)-44*x(453)-47*x(454)-41*x(455)-40*x(456)-48*x(457)-41*x(458)-46*x(459)-49*x(460)-46*x(461)-40*x(462)-42*x(463)-42*x(464)-41*x(465)-41*x(466)-44*x(467)-41*x(468)-43*x(469)-41*x(470)-49*x(471)-48*x(472)-48*x(473)-41*x(474)-46*x(475)-41*x(476)-42*x(477)-42*x(478)-43*x(479)-47*x(480)-43*x(481)-46*x(482)-49*x(483)-41*x(484)-46*x(485)-48*x(486)-42*x(487)-43*x(488)-44*x(489)-43*x(490)-47*x(491)-48*x(492)-42*x(493)-50*x(494)-49*x(495)-41*x(496)-46*x(497)-49*x(498)-42*x(499)-43*x(500)-47*x(501)-43*x(502)-41*x(503)-46*x(504)-42*x(505)-43*x(506)-41*x(507)-50*x(508)-49*x(509)-43*x(510)-43*x(511)-40*x(512)-42*x(513)-41*x(514)-41*x(515)-46*x(516)-49*x(517)-41*x(518)-46*x(519)-49*x(520)-42*x(521)-48*x(522)-43*x(523)-42*x(524)-50*x(525)-49*x(526)-41*x(527)-43*x(528)-44*x(529)-43*x(530)-46*x(531)-48*x(532)-42*x(533)-43*x(534)-41*x(535)-42*x(536)-48*x(537)-43*x(538)-42*x(539)-41*x(540)-50*x(541)-46*x(542)-48*x(543)-48*x(544)-43*x(545)-41*x(546)-46*x(547)-41*x(548)-46*x(549)-49*x(550)-43*x(551)-41*x(552)-43*x(553)-44*x(554)-43*x(555)-41*x(556)-49*x(557)-40*x(558);
end
